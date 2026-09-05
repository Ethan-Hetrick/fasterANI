//! Memory-mapped views over cached reference sketches.

use std::{
    fs, io,
    mem::{align_of, size_of},
    os::fd::AsRawFd,
    path::Path,
    ptr::NonNull,
    sync::Arc,
};

use boomphf::Mphf;

use crate::ani::{
    constants::MinimizerKey,
    model::reference::{
        ContigRecord, ReferenceContigs, ReferenceIndex, ReferenceMinimizer, ReferenceSketch,
        SeedHit,
    },
};

/// Memory-mapped minimal-perfect-hash index and the arrays it addresses.
pub(crate) struct MmapReferenceIndex {
    pub(crate) mphf: Mphf<MinimizerKey>,
    pub(crate) mmap: Arc<MmapFile>,
    pub(crate) key_count: usize,
    pub(crate) hit_count: usize,
    pub(crate) slot_keys_offset: usize,
    pub(crate) hit_offsets_offset: usize,
    pub(crate) hit_counts_offset: usize,
    pub(crate) hit_payloads_offset: usize,
}

impl MmapReferenceIndex {
    pub(crate) fn get(&self, minimizer: &MinimizerKey) -> Option<&[SeedHit]> {
        let slot = self.mphf.try_hash(minimizer)? as usize;
        if slot >= self.key_count {
            return None;
        }

        let slot_keys = self.slot_keys();
        if slot_keys[slot] != *minimizer {
            return None;
        }

        let start = self.hit_offsets()[slot] as usize;
        let count = self.hit_counts()[slot] as usize;
        let end = start.checked_add(count)?;
        self.hit_payloads().get(start..end)
    }

    pub(crate) fn hit_range_by_slot(
        &self,
        slot: usize,
        minimizer: &MinimizerKey,
    ) -> Option<(u32, u32)> {
        if slot >= self.key_count {
            return None;
        }
        if self.slot_keys()[slot] != *minimizer {
            return None;
        }
        Some((self.hit_offsets()[slot], self.hit_counts()[slot]))
    }

    pub(crate) fn hit_payload_range(&self, offset: u32, count: u32) -> Option<&[SeedHit]> {
        let start = offset as usize;
        let count = count as usize;
        let end = start.checked_add(count)?;
        self.hit_payloads().get(start..end)
    }

    pub(crate) fn slot_keys(&self) -> &[MinimizerKey] {
        mmap_slice_at(&self.mmap, self.slot_keys_offset, self.key_count)
    }

    pub(crate) fn hit_offsets(&self) -> &[u32] {
        mmap_slice_at(&self.mmap, self.hit_offsets_offset, self.key_count)
    }

    pub(crate) fn hit_counts(&self) -> &[u32] {
        mmap_slice_at(&self.mmap, self.hit_counts_offset, self.key_count)
    }

    pub(crate) fn hit_payloads(&self) -> &[SeedHit] {
        mmap_slice_at(&self.mmap, self.hit_payloads_offset, self.hit_count)
    }
}

/// Memory-mapped contig records and the flat minimizer array they address.
pub(crate) struct MmapReferenceContigs {
    pub(crate) mmap: Arc<MmapFile>,
    pub(crate) contig_count: usize,
    pub(crate) reference_minimizer_count: usize,
    pub(crate) contig_records_offset: usize,
    pub(crate) reference_minimizers_offset: usize,
}

impl MmapReferenceContigs {
    pub(crate) fn records(&self) -> &[ContigRecord] {
        mmap_slice_at(&self.mmap, self.contig_records_offset, self.contig_count)
    }

    pub(crate) fn reference_minimizers(&self) -> &[ReferenceMinimizer] {
        mmap_slice_at(
            &self.mmap,
            self.reference_minimizers_offset,
            self.reference_minimizer_count,
        )
    }
}

impl ReferenceContigs {
    pub(crate) fn len(&self) -> usize {
        match self {
            Self::Owned(contigs) => contigs.len(),
            Self::Mmap(contigs) => contigs.contig_count,
        }
    }

    pub(crate) fn file_id(&self, contig_id: usize) -> Option<usize> {
        match self {
            Self::Owned(contigs) => contigs.get(contig_id).map(|contig| contig.file_id),
            Self::Mmap(contigs) => contigs
                .records()
                .get(contig_id)
                .map(|record| record.file_id as usize),
        }
    }

    pub(crate) fn minimizers(&self, contig_id: usize) -> Option<&[ReferenceMinimizer]> {
        match self {
            Self::Owned(contigs) => contigs
                .get(contig_id)
                .map(|contig| contig.minimizers.as_slice()),
            Self::Mmap(contigs) => {
                let record: ContigRecord = *contigs.records().get(contig_id)?;
                let start: usize = record.minimizer_offset as usize;
                let count: usize = record.minimizer_count as usize;
                let end: usize = start.checked_add(count)?;
                contigs.reference_minimizers().get(start..end)
            }
        }
    }

    #[cfg(debug_assertions)]
    pub(crate) fn total_minimizers(&self) -> usize {
        match self {
            Self::Owned(contigs) => contigs.iter().map(|contig| contig.minimizers.len()).sum(),
            Self::Mmap(contigs) => contigs.reference_minimizer_count,
        }
    }

    #[cfg(debug_assertions)]
    pub(crate) fn owned_minimizer_capacity_bytes(&self) -> usize {
        match self {
            Self::Owned(contigs) => contigs
                .iter()
                .map(|contig| contig.minimizers.capacity() * size_of::<ReferenceMinimizer>())
                .sum(),
            Self::Mmap(_) => 0,
        }
    }
}

fn mmap_slice_at<T>(mmap: &MmapFile, offset: usize, count: usize) -> &[T] {
    let byte_len: usize = count
        .checked_mul(size_of::<T>())
        .expect("mmap slice length overflow");
    let end: usize = offset
        .checked_add(byte_len)
        .expect("mmap slice offset overflow");
    assert!(end <= mmap.as_slice().len());
    assert_eq!(
        (mmap.as_slice().as_ptr() as usize + offset) % align_of::<T>(),
        0
    );

    unsafe { std::slice::from_raw_parts(mmap.as_slice().as_ptr().add(offset).cast::<T>(), count) }
}

/// Read-only memory map wrapper for cached reference sketches.
pub(crate) struct MmapFile {
    pub(crate) ptr: NonNull<u8>,
    pub(crate) len: usize,
}

unsafe impl Send for MmapFile {}
unsafe impl Sync for MmapFile {}

impl MmapFile {
    pub(crate) fn open(path: &Path) -> io::Result<Self> {
        let file: fs::File = fs::File::open(path)?;
        let len: usize = file.metadata()?.len() as usize;
        if len == 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "cannot mmap an empty reference sketch",
            ));
        }

        let ptr: *mut libc::c_void = unsafe {
            libc::mmap(
                std::ptr::null_mut(),
                len,
                libc::PROT_READ,
                libc::MAP_PRIVATE,
                file.as_raw_fd(),
                0,
            )
        };

        if ptr == libc::MAP_FAILED {
            return Err(io::Error::last_os_error());
        }

        // Disable kernel readahead. Access pattern is random across minimizer
        // positions; readahead fetches pages that will not be used and evicts
        // pages that will be needed.
        unsafe {
            libc::madvise(ptr, len, libc::MADV_RANDOM);
        }

        Ok(Self {
            ptr: NonNull::new(ptr.cast::<u8>()).expect("mmap returned null"),
            len,
        })
    }

    pub(crate) fn as_slice(&self) -> &[u8] {
        unsafe { std::slice::from_raw_parts(self.ptr.as_ptr(), self.len) }
    }

    pub(crate) fn prefetch_sequential(&self) {
        unsafe {
            libc::madvise(
                self.ptr.as_ptr().cast::<libc::c_void>(),
                self.len,
                libc::MADV_SEQUENTIAL,
            );
        }
    }
}

impl ReferenceSketch {
    pub(crate) fn prefetch_sequential(&self) {
        if let ReferenceIndex::Mphf(index) = &self.index {
            index.mmap.prefetch_sequential();
        }
        if let ReferenceContigs::Mmap(contigs) = &self.contigs {
            contigs.mmap.prefetch_sequential();
        }
    }
}

impl Drop for MmapFile {
    fn drop(&mut self) {
        unsafe {
            libc::munmap(self.ptr.as_ptr().cast::<libc::c_void>(), self.len);
        }
    }
}
