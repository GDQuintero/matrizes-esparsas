module gustavo

    !TYPE QUE EMPACOTA UM VETOR ESPARSO
    type PackedVector
        integer :: nz, NFull
        real, allocatable :: Value(:)
        integer, allocatable :: Vector_Index(:)
    end type 
    
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
    
    type Pivot
        integer :: Row, Col
    end type
    
    real :: VectorDensity = 0.25, MatrixDensity = 0.6
    real :: u = 0.25 !PARAMETRO LIMITE (THRESHOLD PIVOTING)
    
    contains
    !================================================================================================
    ! FUNCAO PARA EMPACOTAR UM VETOR ESPARSO 
    !================================================================================================
    function Gather(FullVector,VectorDensity)
        implicit none
        
        integer :: i, j, NonZero, Density, n, VectorDensity
        real :: FullVector(:)
        type(PackedVector) :: Gather
        
        NonZero = 0; n = size(FullVector); j = 0; Density = int(n*VectorDensity)+1
        
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
    function GatherCol(FullMatrix,MatrixDensity)
        implicit none
        
        integer :: i, j, k, Density, m, n, NonZero
        real :: FullMatrix(:,:), MatrixDensity
        type(ColPacked) :: GatherCol
        
        m = size(FullMatrix(:,1)); n = size(FullMatrix(1,:)); Density = int(m*n*MatrixDensity)+1; NonZero = 0; k = 0
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
    function GatherRow(FullMatrix,MatrixDensity)
        implicit none
        
        integer :: i, j, k, Density, m, n, NonZero
        real :: FullMatrix(:,:), MatrixDensity
        type(RowPacked) :: GatherRow
        
        m = size(FullMatrix(:,1)); n = size(FullMatrix(1,:))
        Density = int(m*n*MatrixDensity)+1; NonZero = 0; k = 0
        allocate(GatherRow%Len_Row(n),GatherRow%Row_Start(n))
        allocate(GatherRow%Col_Index(Density),GatherRow%Value(Density))
        
        do i = 1, m
            GatherRow%Row_Start(i) = k + 1
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
    
    !================================================================================================
    ! SOMA DE DUA LINHAS DE UMA MATRIZ EMPACOTADA COMO COLECAO DE LINHAS
    !================================================================================================
    function RowSumColPacked(A,ind1,ind2,alpha,w)
        implicit none
        
        integer :: ind1, ind2, i, j, n, m, zeros, k, NonZero
        real :: alpha, w(:)
        real, allocatable :: aux1(:,:), aux2(:,:)
        type(RowPacked) :: A, RowSumColPacked
        
        m = size(A%Len_Row); n = A%Row_Start(m) + A%Len_Row(m) - A%Row_Start(ind1+1) 
        

        allocate(aux1(2,n))
        w = 0.d0; j = A%Len_Row(ind1) + A%Row_Start(ind1); zeros = 0; k = 0; NonZero = 0
        aux1(1,:) = A%Col_Index(A%Row_Start(ind1+1):A%Row_Start(m)+A%Len_Row(m)-1)
        aux1(2,:) = A%Value(A%Row_Start(ind1+1):A%Row_Start(m)+A%Len_Row(m)-1)
        
        do i = A%Row_Start(ind2), A%Row_Start(ind2) + A%Len_Row(ind2) - 1
            w(A%Col_Index(i)) = A%Value(i)
        enddo
        
        do i = A%Row_Start(ind1), A%Row_Start(ind1) + A%Len_Row(ind1) - 1
            if (w(A%Col_Index(i)) .ne. 0) then
                A%Value(i) = A%Value(i) + alpha*w(A%Col_Index(i))
                w(A%Col_Index(i)) = 0.d0
                
                if (A%Value(i) .eq. 0) then
                    zeros = zeros + 1
                endif
            endif
        enddo
        
        do i = A%Row_Start(ind2), A%Row_Start(ind2) + A%Len_Row(ind2) - 1
            if (W(A%Col_Index(i)) .ne. 0) then
                A%Col_Index(j) = A%Col_Index(i)
                A%Value(j) = alpha*w(A%Col_Index(i))
                w(A%Col_Index(i)) = 0.d0
                j = j + 1; NonZero = NonZero + 1
                if (A%Value(j) .eq. 0) then
                    zeros = zeros + 1
                endif
            endif
        enddo
        
        A%Len_Row(ind1) = A%Len_Row(ind1) + NonZero
        A%Row_Start(ind1+1:) = A%Row_Start(ind1+1:) + NonZero 
       
        if (zeros .ne. 0) then
            allocate(aux2(2,j-zeros-1))
            do i = A%Row_Start(ind1), j-1
                if (A%Value(i) .ne. 0) then
                    k = k + 1
                    aux2(1,k) = A%Col_Index(i)
                    aux2(2,k) = A%Value(i)
                endif
            enddo
            A%Col_Index(A%Row_Start(ind1):) = aux2(1,:)
            A%Value(A%Row_Start(ind1):) = aux2(2,:)
            A%Col_Index(A%Len_Row(ind1)+1:) = aux1(1,:)
            A%Value(A%Len_Row(ind1)+1:) = aux1(2,:)
        endif
        
        A%Col_Index(A%Row_Start(ind1+1):) = aux1(1,:)
        A%Value(A%Row_Start(ind1+1):) = aux1(2,:)
        RowSumColPacked = A
    end function RowSumColPacked
    
    !================================================================================================
    ! CRITERIO DE GRADO MINIMO
    !================================================================================================
    function MinDeg(A)
        implicit none
        
        integer :: i, j, k, l, rk, n
        integer, allocatable :: ind(:)
        type(RowPacked) :: A
        type(Pivot) :: MinDeg  
        real :: akj, akk
        
        rk = minval(A%Len_Row); n = size(A%Len_Row); j = 0; k = 0; l = 0
        
        do i = 1, n
            if (A%Len_Row(i) .eq. rk) then
                l = l + 1
            endif
        enddo
        
        allocate(ind(l))
        
        do i = 1, n
            j = j + 1
            if (A%Len_Row(i) .eq. rk) then
                k = k + 1
                ind(k) = j
            endif 
        enddo
        
        if (j .eq. 1) then
            !MinDeg%Value = A%Value(A%Row_Start(ind(1)))
            MinDeg%Row = ind(1)
            MinDeg%Col = MinDeg%Row
            return
        endif
        
        k = 0; j = 0
        
        do i = 1, l
            do k = 1, rk
                if (A%Col_Index(A%Row_Start(ind(i)) + k - 1) .eq. ind(i)) then
                    akk = A%Value(A%Row_Start(ind(i)) + k - 1)
                    exit
                endif
            enddo
            do j = 1, rk
                akj = A%Value(A%Row_Start(ind(i)) + j - 1)
                if (A%Col_Index(A%Row_Start(ind(i)) + j - 1) .ne. ind(i) .and. abs(akk) .ge. u*abs(akj)) then
                    !MinDeg%Value = akk
                    MinDeg%Row = ind(i)
                    MinDeg%Col = MinDeg%Row
                    return
                endif
            enddo
        enddo
        
    end function MinDeg
    
end module
