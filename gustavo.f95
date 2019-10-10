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
        real :: Value
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
    function GatherCol(FullMatrix)
        implicit none
        
        integer :: i, j, k, Density, m, n, NonZero
        real :: FullMatrix(:,:)
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
    function GatherRow(FullMatrix)
        implicit none
        
        integer :: i, j, k, Density, m, n, NonZero
        real :: FullMatrix(:,:)
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
    ! DESEMPACOTAMENTO
    !================================================================================================    
    function Unpaking(X)
        implicit none
        
        integer :: m, n, NonZero = 0, i, j
        real, allocatable :: Unpaking(:,:)
        type(RowPacked) :: X
        
        n = size(X%Len_Row); m = maxval(X%Col_Index)
        allocate(Unpaking(m,n))
        
        do i = 1, n
            do j = X%Row_Start(i), X%Row_Start(i) + X%Len_Row(i) - 1
                Unpaking(i,X%Col_Index(j)) = X%Value(j)  
            enddo
        enddo
                
    end function Unpaking
    
    !================================================================================================
    ! EMPACOTAMENTO - FORMA LINHA A FORMA COLUNA
    !================================================================================================ 
    function PackRowCol(A)
        implicit none
        
        type(RowPacked) :: A
        type(ColPacked) :: PackRowCol
        integer :: n=0, m=0, k=0, i=0, j=0, l=0
        
        m = size(A%Len_Row)
        allocate(PackRowCol%Len_Col(m),PackRowCol%Col_Start(m))
        PackRowCol%Len_Col = 0; PackRowCol%Col_Start = 1
                
        do i = 1, m
            k = k + A%Len_Row(i)
        enddo
        
        allocate(PackRowCol%Row_Index(k),PackRowCol%Value(k))
        PackRowCol%Row_Index = 0; PackRowCol%Value = 0
        
        do i = 1, k
            PackRowCol%Len_Col(A%Col_Index(i)) = PackRowCol%Len_Col(A%Col_Index(i)) + 1
        enddo
        
        do i = 2, m
            PackRowCol%Col_Start(i) = PackRowCol%Col_Start(i-1) + PackRowCol%Len_Col(i-1)
        enddo
        k = 0
        do i = 1, m
            do j = 1, A%Len_Row(i)
                l = l + 1
                do while (PackRowCol%Row_Index(PackRowCol%Col_Start(A%Col_Index(l))+k) .ne. 0)
                    k = k + 1
                enddo
                PackRowCol%Row_Index(PackRowCol%Col_Start(A%Col_Index(l))+k) = i
                PackRowCol%Value(PackRowCol%Col_Start(A%Col_Index(l))+k) = A%Value(l)
                k = 0
            enddo
        enddo
        
    end function PackRowCol
    
    !================================================================================================
    ! EMPACOTAMENTO - FORMA COLUNA A FORMA LINHA
    !================================================================================================ 
    function PackColRow(A)
        implicit none
        
        type(RowPacked) :: PackColRow
        type(ColPacked) :: A
        integer :: n=0, m=0, k=0, i=0, j=0, l=0
        
        m = size(A%Len_Col)
        allocate(PackColRow%Len_Row(m),PackColRow%Row_Start(m))
        PackColRow%Len_Row = 0; PackColRow%Row_Start = 1
                
        do i = 1, m
            k = k + A%Len_Col(i)
        enddo
        
        allocate(PackColRow%Col_Index(k),PackColRow%Value(k))
        PackColRow%Col_Index = 0; PackColRow%Value = 0
        
        do i = 1, k
            PackColRow%Len_Row(A%Row_Index(i)) = PackColRow%Len_Row(A%Row_Index(i)) + 1
        enddo
        
        do i = 2, m
            PackColRow%Row_Start(i) = PackColRow%Row_Start(i-1) + PackColRow%Len_Row(i-1)
        enddo
        k = 0
        do i = 1, m
            do j = 1, A%Len_Col(i)
                l = l + 1
                do while (PackColRow%Col_Index(PackColRow%Row_Start(A%Row_Index(l))+k) .ne. 0)
                    k = k + 1
                enddo
                PackColRow%Col_Index(PackColRow%Row_Start(A%Row_Index(l))+k) = i
                PackColRow%Value(PackColRow%Row_Start(A%Row_Index(l))+k) = A%Value(l)
                k = 0
            enddo
        enddo
        
    end function PackColRow
    
    !================================================================================================
    ! SOMA DE DUAS LINHAS DE UMA MATRIZ EMPACOTADA COMO COLECAO DE LINHAS (x + alpha*y)
    !================================================================================================
    subroutine RRowSumColPacked(A,ind1,ind2,alpha,w)
        implicit none
        
        integer :: ind1, ind2, i, j, n, m, zeros, k, NonZero
        real :: alpha, w(:)
        real, allocatable :: aux1(:,:), aux2(:,:)
        type(RowPacked) :: A
        
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
        
        A%Len_Row(ind1) = A%Len_Row(ind1) + NonZero - zeros
        A%Row_Start(ind1+1:) = A%Row_Start(ind1+1:) + NonZero
       
        if (zeros .ne. 0) then
            allocate(aux2(2,A%Len_Row(ind1)+NonZero-zeros))
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
    end subroutine RRowSumColPacked
    
    !================================================================================================
    ! SOMA DE DUAS LINHAS DE UMA MATRIZ EMPACOTADA COMO COLECAO DE LINHAS (alpha*x + y)
    !================================================================================================
    
    subroutine RowSumRowPacked(A,ind1,ind2,alpha,w)
        implicit none
        
        integer :: ind1, ind2, i, j, n, m, zeros, k, NonZero
        real :: alpha, w(:)
        real, allocatable :: aux1(:,:), aux2(:,:)
        type(RowPacked) :: A
        
        m = size(A%Len_Row); n = A%Row_Start(m) + A%Len_Row(m) - A%Row_Start(ind1+1) 
        
        allocate(aux1(2,n))
        j = A%Len_Row(ind1) + A%Row_Start(ind1); zeros = 0; k = 0; NonZero = 0
        aux1(1,:) = A%Col_Index(A%Row_Start(ind1+1):A%Row_Start(m)+A%Len_Row(m)-1)
        aux1(2,:) = A%Value(A%Row_Start(ind1+1):A%Row_Start(m)+A%Len_Row(m)-1)
        
        do i = A%Row_Start(ind2), A%Row_Start(ind2) + A%Len_Row(ind2) - 1
            w(A%Col_Index(i)) = A%Value(i)
        enddo

        do i = A%Row_Start(ind1), A%Row_Start(ind1) + A%Len_Row(ind1) - 1
            A%Value(i) = alpha*A%Value(i)
        enddo
        
        do i = A%Row_Start(ind1), A%Row_Start(ind1) + A%Len_Row(ind1) - 1
            if (w(A%Col_Index(i)) .ne. 0) then
                A%Value(i) = A%Value(i) + w(A%Col_Index(i))
                w(A%Col_Index(i)) = 0.d0
                if (A%Value(i) .eq. 0) then
                    zeros = zeros + 1
                endif
            endif
        enddo
        
        do i = A%Row_Start(ind2), A%Row_Start(ind2) + A%Len_Row(ind2) - 1
            if (W(A%Col_Index(i)) .ne. 0) then
                A%Col_Index(j) = A%Col_Index(i)
                A%Value(j) = w(A%Col_Index(i))
                w(A%Col_Index(i)) = 0.d0
                j = j + 1; NonZero = NonZero + 1
                if (A%Value(j) .eq. 0) then
                    zeros = zeros + 1
                endif
            endif
        enddo
        
        A%Len_Row(ind1) = A%Len_Row(ind1) + NonZero - zeros
        A%Row_Start(ind1+1:) = A%Row_Start(ind1+1:) + NonZero - zeros
       
        if (zeros .ne. 0) then
            allocate(aux2(2,A%Len_Row(ind1)))
            do i = A%Row_Start(ind1), j-1
                if (A%Value(i) .ne. 0) then
                    k = k + 1
                    aux2(1,k) = A%Col_Index(i)
                    aux2(2,k) = A%Value(i)
                endif
            enddo
            A%Col_Index(A%Row_Start(ind1):) = aux2(1,:)
            A%Value(A%Row_Start(ind1):) = aux2(2,:)
            A%Col_Index(A%Row_Start(ind1+1):) = aux1(1,:)
            A%Value(A%Row_Start(ind1+1):) = aux1(2,:)
        endif
        A%Col_Index(A%Row_Start(ind1+1):) = aux1(1,:)
        A%Value(A%Row_Start(ind1+1):) = aux1(2,:)
    end subroutine RowSumRowPacked
    
    !================================================================================================
    ! CRITERIO DE GRADO MINIMO
    !================================================================================================
    function MinDeg(A)
        implicit none
        
        integer :: i = 0, r_i = 0
        type(RowPacked) :: A
        type(Pivot) :: MinDeg
        
        !Contamos los elementos no nulos de la primera fila
        r_i = A%Len_Row(1)
        MinDeg%Row = 1
        MinDeg%Col = 1
        
        !Buscamos min r_i, y seleccionmos el elemento diagonal como pivote
        do i = 2, size(A%Len_Row)
            if (A%Len_Row(i) .lt. r_i) then
                r_i = A%Len_Row(i)
                MinDeg%Row = i
                MinDeg%Col = i
            endif
        enddo
        
        !Buscamos el valor del pivote
        do i = 1, r_i
            if (A%Col_Index(A%Row_Start(MinDeg%Row)+i-1) .eq. MinDeg%Row) then
                MinDeg%Value = A%Value(A%Row_Start(MinDeg%Row)+i-1)
            endif
        enddo
    end function MinDeg
    
end module
