Module Tools
    implicit none
    
    !TYPE QUE EMPACOTA UM VETOR ESPARSO
    type PackedVector
        integer :: NZ, NFull
        real, allocatable :: Value(:)
        integer, allocatable :: VectorIndex(:)
    end type PackedVector
    
    !TYPE QUE EMPACOTA UMA MATRIZ NO ESQUEMA COORDENADAS
    type PackedMatrixEntry
        integer :: NNZ
        integer, allocatable :: IRN(:), JCN(:)
        real, allocatable :: Value(:)
    end type PackedMatrixEntry
    
    !TYPE QUE EMPACOTA UMA MATRIZ COMO COLECAO DE COLUNAS
    type PackedMatrixCol
        integer, allocatable :: LenCol(:), ColStart(:), RowIndex(:)
        real, allocatable :: Value(:)
    end type PackedMatrixCol
    
    !TYPE QUE EMPACOTA UMA MATRIZ COMO COLECAO DE LINHAS
    type PackedMatrixRow
        integer, allocatable :: LenRow(:), RowStart(:), ColIndex(:)
        real, allocatable :: Value(:)
    end type PackedMatrixRow
    
    
end Module Tools
