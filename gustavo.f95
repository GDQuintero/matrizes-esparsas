module Gustavo
    !TYPE QUE EMPACOTA UM VETOR ESPARSO
    type PackedVector
        integer :: nz, NFull
        real, allocatable :: Value(:)
        integer, allocatable :: Vector_Index(:)
    end type PackedVector
    
    !TYPE QUE EMPACOTA UMA MATRIZ NO ESQUEMA COORDENADAS
    type EntryPacked
        integer :: NNZ
        integer, allocatable :: Row_Index(:), Col_Index(:)
        real, allocatable :: Value(:)
    end type 
    
    !TYPE QUE EMPACOTA UMA MATRIZ COMO COLECAO DE COLUNAS
    type ColPacked
        integer, allocatable :: Len_Col(:), Col_Start(:), Row_Index(:)
        real, allocatable :: Value(:)
    end type 
    
    !TYPE QUE EMPACOTA UMA MATRIZ COMO COLECAO DE LINHAS
    type RowPacked
        integer, allocatable :: Len_Row(:), Row_Start(:), Col_Index(:)
        real, allocatable :: Value(:)
    end type 
    
    contains
    !================================================================================================
    ! FUNCAO PARA EMPACOTAR UM VETOR ESPARSO 
    !================================================================================================
    function Gather(FullVector)
        implicit none
        
        integer :: i, j, NonZero, Density, n
        real :: FullVector(:)
        type(PackedVector) :: Gather
        
        NonZero = 0; n = size(FullVector); j = 0; Density = int(n*0.5d0)+1
        
        do i = 1, n
            if (FullVector(i) .ne. 0.d0) then
                NonZero = NonZero + 1
            endif
        enddo
        
        Gather%NZ = NonZero; Gather%Nfull = n
        allocate(Gather%Vector_Index(Density),Gather%Value(Density))
        
        do i = 1, n
            if (FullVector(i) .ne. 0.d0) then
                j = j + 1
                Gather%Vector_Index(j) = i
                Gather%Value(j) = FullVector(i)
            endif 
        enddo
    end function Gather
end module
