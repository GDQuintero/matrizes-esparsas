module Tools
    implicit none
    
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
    
!     !================================================================================================
!     ! FUNCAO PARA EMPACOTAR UM VETOR ESPARSO 
!     !================================================================================================
!     function Gather(FullVector)
!         implicit none
!         
!         integer :: i, j, NonZero, Density, n
!         real :: FullVector(:)
!         type(PackedVector) :: Gather
!         
!         NonZero = 0; n = size(FullVector); j = 0; Density = int(n*0.5d0)+1
!         
!         do i = 1, n
!             if (FullVector(i) .ne. 0.d0) then
!                 NonZero = NonZero + 1
!             endif
!         enddo
!         
!         Gather%NZ = NonZero; Gather%Nfull = n
!         allocate(Gather%Vector_Index(Density),Gather%Value(Density))
!         
!         do i = 1, n
!             if (FullVector(i) .ne. 0.d0) then
!                 j = j + 1
!                 Gather%Vector_Index(j) = i
!                 Gather%Value(j) = FullVector(i)
!             endif 
!         enddo
!     end function Gather
!     
!     !================================================================================================
!     ! FUNCAO PARA EMPACOTAR UMA MATRIZ ESPARSA (COORD. SCHEME)
!     !================================================================================================
!     function GatherMatrixEntry(FullMatrix)
!         implicit none 
!         
!         integer :: i, j, k, NonZero, Density, n, m
!         real :: FullMatrix(:,:)
!         type(PackedMatrixEntry) :: GatherMatrixEntry
!         
!         NonZero = 0; k = 0
!         m = size(FullMatrix(:,1)); n = size(FullMatrix(1,:)); Density = int(m*n*0.5)+1
!         
!         do i = 1, n
!             do j = 1, m
!                 if (FullMatrix(i,j) .ne. 0) then
!                     NonZero = NonZero + 1
!                 endif
!             enddo
!         enddo
! 
!         allocate(GatherMatrixEntry%Row_Index(Density),GatherMatrixEntry%Col_Index(Density),GatherMatrixEntry%Value(Density))
!         GatherMatrixEntry%NNZ = NonZero
!         
!         do i = 1, n
!             do j = 1, m
!                 if (FullMatrix(i,j) .ne. 0) then
!                     k = k + 1
!                     GatherMatrixEntry%Row_Index(k) = i
!                     GatherMatrixEntry%Col_Index(k) = j
!                     GatherMatrixEntry%Value(k) = FullMatrix(i,j)
!                 endif
!             enddo
!         enddo
!     end function GatherMatrixEntry
!     
!     !================================================================================================
!     ! FUNCAO PARA EMPACOTAR UMA MATRIZ ESPARSA (POR COLUNAS)
!     !================================================================================================
!     function GatherMatrixCol(FullMatrix)
!         implicit none
!         
!         integer :: i, j, k, Density, m, n, NonZero
!         real :: FullMatrix(:,:)
!         type(PackedMatrixCol) :: GatherMatrixCol
!         
!         m = size(FullMatrix(:,1)); n = size(FullMatrix(1,:)); Density = int(m*n*0.5)+1; NonZero = 0; k = 0
!         allocate(GatherMatrixCol%Len_Col(n),GatherMatrixCol%Col_Start(n))
!         allocate(GatherMatrixCol%Row_Index(Density),GatherMatrixCol%Value(Density))
!         
!         do j = 1, n
!             GatherMatrixCol%Col_Start(j) = 1 + k
!             do i = 1,m
!                 if (FullMatrix(i,j) .ne. 0) then
!                     k = k + 1
!                     GatherMatrixCol%Row_Index(k) = i
!                     GatherMatrixCol%Value(k) = FullMatrix(i,j)
!                     NonZero = NonZero + 1
!                 endif
!             enddo
!             GatherMatrixCol%Len_Col(j) = NonZero
!             NonZero = 0
!         enddo
!     end function GatherMatrixCol
!     
!     !================================================================================================
!     ! FUNCAO PARA EMPACOTAR UMA MATRIZ ESPARSA (POR LINHAS)
!     !================================================================================================
!     function GatherMatrixRow(FullMatrix)
!         implicit none
!         
!         integer :: i, j, k, Density, m, n, NonZero
!         real :: FullMatrix(:,:)
!         type(PackedMatrixRow) :: GatherMatrixRow
!         
!         m = size(FullMatrix(:,1)); n = size(FullMatrix(1,:)); Density = int(m*n*0.5)+1; NonZero = 0; k = 0
!         allocate(GatherMatrixRow%Len_Row(n),GatherMatrixRow%Row_Start(n))
!         allocate(GatherMatrixRow%Col_Index(Density),GatherMatrixRow%Value(Density))
!         
!         do i = 1, m
!             GatherMatrixRow%Row_Start(i) = 1 + k
!             do j = 1,n
!                 if (FullMatrix(i,j) .ne. 0) then
!                     k = k + 1
!                     GatherMatrixRow%Col_Index(k) = j 
!                     GatherMatrixRow%Value(k) = FullMatrix(i,j)
!                     NonZero = NonZero + 1
!                 endif
!             enddo
!             GatherMatrixRow%Len_Row(i) = NonZero
!             NonZero = 0
!         enddo
!     end function GatherMatrixRow
    
end module Tools
