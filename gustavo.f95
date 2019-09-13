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
    
    !================================================================================================
    ! FUNCAO PARA EMPACOTAR UMA MATRIZ ESPARSA (POR COLUNAS)
    !================================================================================================
    function GatherCol(FullMatrix)
        implicit none
        
        integer :: i, j, k, Density, m, n, NonZero
        real :: FullMatrix(:,:)
        type(ColPacked) :: GatherCol
        
        m = size(FullMatrix(:,1)); n = size(FullMatrix(1,:)); Density = int(m*n*0.5)+1; NonZero = 0; k = 0
        allocate(GatherCol%Len_Col(n),GatherCol%Col_Start(n))
        allocate(GatherCol%Row_Index(Density),GatherCol%Value(Density))
        
        do j = 1, n
            GatherCol%Col_Start(j) = 1 + k
            do i = 1,m
                if (FullMatrix(i,j) .ne. 0) then
                    k = k + 1
                    GatherCol%Row_Index(k) = i
                    GatherCol%Value(k) = FullMatrix(i,j)
                    NonZero = NonZero + 1
                endif
            enddo
            GatherCol%Len_Col(j) = NonZero
            NonZero = 0
        enddo
    end function GatherCol
    
    !================================================================================================
    ! FUNCAO PARA EMPACOTAR UMA MATRIZ ESPARSA (POR LINHAS)
    !================================================================================================
    function GatherRow(FullMatrix)
        implicit none
        
        integer :: i, j, k, Density, m, n, NonZero
        real :: FullMatrix(:,:)
        type(RowPacked) :: GatherRow
        
        m = size(FullMatrix(:,1)); n = size(FullMatrix(1,:)); Density = int(m*n*0.5)+1; NonZero = 0; k = 0
        allocate(GatherRow%Len_Row(n),GatherRow%Row_Start(n))
        allocate(GatherRow%Col_Index(Density),GatherRow%Value(Density))
        
        do i = 1, m
            GatherRow%Row_Start(i) = 1 + k
            do j = 1,n
                if (FullMatrix(i,j) .ne. 0) then
                    k = k + 1
                    GatherRow%Col_Index(k) = j 
                    GatherRow%Value(k) = FullMatrix(i,j)
                    NonZero = NonZero + 1
                endif
            enddo
            GatherRow%Len_Row(i) = NonZero
            NonZero = 0
        enddo
    end function GatherRow
    
    
    
end module
