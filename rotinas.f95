Module Tools
    implicit none
    
    type PackedVector
        integer :: NZ, NFull
        real, allocatable :: Value(:)
        integer, allocatable :: VectorIndex(:)
    end type PackedVector

    type PackedMatrixEntry
        integer :: NNZ
        integer, allocatable :: IRN(:), JCN(:)
        real, allocatable :: Value(:)
    end type PackedMatrixEntry
    
    type PackedMatrixCol
        integer, allocatable :: LenCol(:), ColStart(:), RowIndex(:)
        real, allocatable :: Value(:)
    end type PackedMatrixCol

    type PackedMatrixRow
        integer, allocatable :: LenRow(:), RowStart(:), ColIndex(:)
        real, allocatable :: Value(:)
    end type PackedMatrixRow
    
    
end Module Tools
