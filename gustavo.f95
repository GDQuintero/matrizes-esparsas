module gustavo

    !TYPE QUE EMPACOTA UM VETOR ESPARSO
    type PackedVector
        integer :: nz, NFull
        real, allocatable :: Value(:)
        integer, allocatable :: Vector_Index(:)
    end type 
    
    !TYPE QUE EMPACOTA UMA MATRIZ NO ESQUEMA COORDENADAS
    type CoordPacked
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
        integer :: n, m
        integer, allocatable :: Len_Row(:), Row_Start(:), Col_Index(:)
        real, allocatable :: Value(:)
    end type 
    
    !TYPE PIVO PARA GRAU MINIMO
    type PivotMD
        integer :: Row, Col, rmin
        real :: Value
    end type
    
    !TYPE PIVO
    type Pivot
        integer :: Row, Col
        real :: Value
    end type

    !ESTRUTURA TIPO GUSTAVSON
    type GustavsonPacked
        integer, allocatable :: Len_Row(:), Len_Col(:), Row_Start(:), Col_Start(:), Row_Index(:), Col_Index(:)
        real, allocatable :: Value(:)
    end type
    
    real :: VectorDensity = 0.25, MatrixDensity = 0.5
    real :: u = 0.25 !PARAMETRO LIMITE (THRESHOLD PIVOTING)
    
    contains

    !================================================================================================
    ! FUNCAO PARA EMPACOTAR UMA MATRIZ ESPARSA (FORMATO GUSTAVSON)
    !================================================================================================
    function GatherGustavson(FullMatrix)
        implicit none
        
        integer :: i, j, k, Density, m, n, NonZero
        real :: FullMatrix(:,:)
        type(GustavsonPacked) :: GatherGustavson
        
        m = size(FullMatrix(:,1)); n = size(FullMatrix(1,:))
        Density = int(m*n*MatrixDensity)+1; NonZero = 0; k = 0
        allocate(GatherGustavson%Len_Row(n),GatherGustavson%Len_Col(m))
        allocate(GatherGustavson%Row_Start(n),GatherGustavson%Col_Start(m))
        allocate(GatherGustavson%Row_Index(Density),GatherGustavson%Col_Index(Density))
        allocate(GatherGustavson%Value(Density))
        
        do i = 1, m
            GatherGustavson%Row_Start(i) = k + 1
            do j = 1,n
                if (FullMatrix(i,j) .ne. 0) then
                    k = k + 1
                    GatherGustavson%Col_Index(k) = j 
                    GatherGustavson%Value(k) = FullMatrix(i,j)
                    NonZero = NonZero + 1
                endif
            enddo
            GatherGustavson%Len_Row(i) = NonZero
            NonZero = 0;
        enddo
        k = 0
        do j = 1, n
            GatherGustavson%Col_Start(j) = k + 1
            do i = 1,m
                if (FullMatrix(i,j) .ne. 0) then
                    k = k + 1
                    GatherGustavson%Row_Index(k) = i
                    NonZero = NonZero + 1
                endif
            enddo
            GatherGustavson%Len_Col(j) = NonZero
            NonZero = 0;
        enddo
        
    end function GatherGustavson

    !================================================================================================
    ! DESEMPACOTAMENTO (FORMATO GUSTAVSON)
    !================================================================================================
    function UnpakinGus(X)
        implicit none
        
        integer :: m, n, NonZero = 0, i, j
        real, allocatable :: UnpakinGus(:,:)
        type(GustavsonPacked) :: X
        
        n = size(X%Len_Row); m = size(X%Len_Col)
        allocate(UnpakinGus(m,n))
        
        do i = 1, n
            do j = X%Row_Start(i), X%Row_Start(i) + X%Len_Row(i) - 1
                UnpakinGus(i,X%Col_Index(j)) = X%Value(j)  
            enddo
        enddo
    end function UnpakinGus

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
        
        m = size(FullMatrix(:,1)); n = size(FullMatrix(1,:))
        Density = int(m*n*MatrixDensity)+1; NonZero = 0; k = 0
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
        GatherRow%n = n; GatherRow%m = m
        Density = int(m*n*MatrixDensity)+1; NonZero = 0; k = 0
        allocate(GatherRow%Len_Row(n),GatherRow%Row_Start(n))
        allocate(GatherRow%Col_Index(Density),GatherRow%Value(Density))
        
        do i = 1, m
            GatherRow%Row_Start(i) = k + 1
            do j = 1, n
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
    ! DESEMPACOTAMENTO (COLECAO DE LINHAS)
    !================================================================================================    
    function Unpaking(X)
        implicit none
        
        integer :: m, n, NonZero = 0, i, j
        real, allocatable :: Unpaking(:,:)
        type(RowPacked) :: X
        
        n = X%n
        allocate(Unpaking(n,n))
        Unpaking = 0.d0
                
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
    subroutine RRowSumRowPacked(A,indx,indy,alpha,w)
        implicit none
        
        integer :: indx, indy, i, j, n, m, zeros, k, NonZero
        real :: alpha, w(:)
        integer, allocatable :: aux11(:), aux21(:)
        real, allocatable :: aux12(:), aux22(:)
        type(RowPacked) :: A
        
        m = size(A%Len_Row); n = A%Row_Start(m) + A%Len_Row(m) - A%Row_Start(indx+1) 
        
        if (indx .ne. m) then!Caso donde el resultado de la suma no es almacenado en la ultima fila
            allocate(aux11(n),aux12(n))
            w = 0.d0; j = A%Len_Row(indx) + A%Row_Start(indx); zeros = 0; k = 0; NonZero = 0
            
            !Guardamos las filas abajo de la fila que se modificara
            aux11 = A%Col_Index(A%Row_Start(indx + 1):A%Row_Start(m) + A%Len_Row(m)-1)
            aux12 = A%Value(A%Row_Start(indx + 1):A%Row_Start(m) + A%Len_Row(m)-1)
            
            !Desempaquetamos una de las filas
            do i = A%Row_Start(indy), A%Row_Start(indy) + A%Len_Row(indy) - 1
                w(A%Col_Index(i)) = A%Value(i)
            enddo
            
            !Sumamos elementos que comparten columnas
            do i = A%Row_Start(indx), A%Row_Start(indx) + A%Len_Row(indx) - 1
                if (w(A%Col_Index(i)) .ne. 0) then
                    A%Value(i) = A%Value(i) + alpha*w(A%Col_Index(i))
                    w(A%Col_Index(i)) = 0.d0
                    
                    if (A%Value(i) .eq. 0) then
                        zeros = zeros + 1
                    endif
                endif
            enddo
            
            !Sumamos elementos restantes
            do i = A%Row_Start(indy), A%Row_Start(indy) + A%Len_Row(indy) - 1
                if (w(A%Col_Index(i)) .ne. 0) then
                    A%Col_Index(j) = A%Col_Index(i)
                    A%Value(j) = alpha*w(A%Col_Index(i))
                    w(A%Col_Index(i)) = 0.d0
                       
                    if (A%Value(j) .eq. 0) then
                        zeros = zeros + 1
                    endif
                    j = j + 1; NonZero = NonZero + 1
                endif
            enddo
            
            !Empujamos a la izquierda si se anulan elementos y a la derecha si surgen nuevas entradas
            if (zeros .ne. 0) then
                allocate(aux21(A%Len_Row(indx) + NonZero - zeros),aux22(A%Len_Row(indx) + NonZero - zeros))
                do i = A%Row_Start(indx), j - 1
                    if (A%Value(i) .ne. 0) then
                        k = k + 1
                        aux21(k) = A%Col_Index(i)
                        aux22(k) = A%Value(i)
                    endif
                enddo
                
                A%Len_Row(indx) = A%Len_Row(indx) + NonZero - zeros
                A%Row_Start(indx + 1:) = A%Row_Start(indx + 1:) + NonZero - zeros
                
                A%Col_Index(A%Row_Start(indx):) = 0
                A%Value(A%Row_Start(indx):) = 0.d0
                
                A%Col_Index(A%Row_Start(indx):) = aux21
                A%Value(A%Row_Start(indx):) = aux22
                
                A%Col_Index(A%Row_Start(indx+1):) = aux11
                A%Value(A%Row_Start(indx+1):) = aux12
                
                return
            else
                A%Len_Row(indx) = A%Len_Row(indx) + NonZero - zeros
                A%Row_Start(indx + 1:) = A%Row_Start(indx + 1:) + NonZero - zeros
                
                A%Col_Index(A%Row_Start(indx + 1):) = aux11
                A%Value(A%Row_Start(indx + 1):) = aux12
            endif
            
        else!Caso donde el resultado de la suma es almacenado en la ultima fila
            
            w = 0.d0; j = A%Len_Row(indx) + A%Row_Start(indx); zeros = 0; k = 0; NonZero = 0
            
            !Desempaquetamos una de las filas
            do i = A%Row_Start(indy), A%Row_Start(indy) + A%Len_Row(indy) - 1
                w(A%Col_Index(i)) = A%Value(i)
            enddo
            
            !Sumamos elementos que comparten columnas
            do i = A%Row_Start(indx), A%Row_Start(indx) + A%Len_Row(indx) - 1
                if (w(A%Col_Index(i)) .ne. 0) then
                    A%Value(i) = A%Value(i) + alpha*w(A%Col_Index(i))
                    w(A%Col_Index(i)) = 0.d0
                    
                    if (A%Value(i) .eq. 0) then
                        zeros = zeros + 1
                    endif
                endif
            enddo
            
            !Sumamos elementos restantes
            do i = A%Row_Start(indy), A%Row_Start(indy) + A%Len_Row(indy) - 1
                if (w(A%Col_Index(i)) .ne. 0) then
                    A%Col_Index(j) = A%Col_Index(i)
                    A%Value(j) = alpha*w(A%Col_Index(i))
                    w(A%Col_Index(i)) = 0.d0
                    
                    if (A%Value(j) .eq. 0) then
                        zeros = zeros + 1
                    endif
                    j = j + 1; NonZero = NonZero + 1
                endif
            enddo
            
            !En este caso solo empujamos a la izquierda si se anulan elementos
            if (zeros .ne. 0) then
                allocate(aux21(A%Len_Row(indx) + NonZero - zeros),aux22(A%Len_Row(indx) + NonZero - zeros))
                do i = A%Row_Start(indx), j - 1
                    if (A%Value(i) .ne. 0) then
                        k = k + 1
                        aux21(k) = A%Col_Index(i)
                        aux22(k) = A%Value(i)
                    endif
                enddo
            
                A%Len_Row(indx) = A%Len_Row(indx) + NonZero - zeros
                A%Row_Start(indx + 1:) = A%Row_Start(indx + 1:) + NonZero - zeros
                
                A%Col_Index(A%Row_Start(indx):) = 0
                A%Value(A%Row_Start(indx):) = 0.d0
                
                A%Col_Index(A%Row_Start(indx):) = aux21
                A%Value(A%Row_Start(indx):) = aux22
                
                return
            else
                A%Len_Row(indx) = A%Len_Row(indx) + NonZero - zeros
                A%Row_Start(indx + 1:) = A%Row_Start(indx + 1:) + NonZero - zeros
            endif
        endif
        
    end subroutine RRowSumRowPacked
    
    !================================================================================================
    ! PRODUTO INTERNO DE DOIS VETORES DENSOS
    !================================================================================================ 
    function ProdVect(x,y)
        implicit none
        
        integer :: i, n
        real :: x(:), y(:), ProdVect
        
        n = size(x) 
        ProdVect = 0.d0
        
        do i = 1, n
            ProdVect = ProdVect + x(i)*y(i)
        enddo
        
    end function ProdVect
    
    !================================================================================================
    ! PRODUTO DE DUAS MATRIZES DENSAS
    !================================================================================================   
    function ProdMat(A,B)
        implicit none
        
        integer :: i, j, n
        real :: A(:,:), B(:,:)
        real, allocatable :: ProdMat(:,:)
        
        n = size(A(1,:)); allocate(ProdMat(n,n))
        
        do i = 1, n
            do j = 1, n
                ProdMat(i,j) = ProdVect(A(i,:),B(:,j))
            enddo
        enddo

    end function ProdMat
    
    !================================================================================================
    ! FUNCAO PARA CRIAR UM VETOR CANONICO
    !================================================================================================ 
    function Canon(ind,n)
        implicit none
    
        integer :: ind, n
        real, allocatable :: Canon(:)
        
        allocate(Canon(n)); Canon = 0.d0; Canon(ind) = 1.d0
        
    end function Canon
        
    !================================================================================================
    ! FUNCAO PARA CRIAR A MATRIZ DE PERMUTACAO
    !================================================================================================ 
    function PerMat(P)
    implicit none
    
    integer :: i, j, n, P(:)
    real, allocatable :: PerMat(:,:)
    
    n = size(P); allocate(PerMat(n,n))
    
    do i = 1, n
        PerMat(i,:) = Canon(P(i),n)
    enddo
    
    end function PerMat
    
    !================================================================================================
    ! FUNCAO PARA IMPRIMIR O PADRAO DE PREENCHIMENTO DE UMA MATRIZ
    !================================================================================================ 
    function Pattern(A,n)
        implicit none
        
        real :: A(:,:)
        character(len=3), allocatable :: Pattern(:,:)
        integer :: i, j, n
        
        allocate(Pattern(n,n))
        Pattern = "   "
        do i = 1, n
            do j = 1, n
                if(abs(A(i,j)) .ge. 10e-4) then
                    Pattern(i,j) = " x "
                else
                    Pattern(i,j) = "   "
                endif
            enddo
        enddo
        
        do i = 1, n
            print*, "|", Pattern(i,:), "|"
        enddo
    end function Pattern
    
    !================================================================================================
    ! FUNCAO PARA IMPRIMIR UMA MATRIZ
    !================================================================================================ 
    subroutine PrintMat(A,n)
        implicit none
        
        real :: A(:,:)
        integer :: i, n
        
        do i = 1, n
            print*, A(i,:)
        enddo
    end subroutine PrintMat
    
    !================================================================================================
    ! FUNCAO PARA EXPORTAR UMA MATRIZ EMPACOTADA EM UM ARQUIVO .DAT
    !================================================================================================ 
    subroutine Export(A,n)
        implicit none
        
        real :: A(:,:)
        integer :: i, n, j
        
        Open(Unit = 10, File = "dados.dat", ACCESS = "SEQUENTIAL")
    
        do i = 1, n
            do j = 1, n
                write(10,*) i, j, A(i,j)
            enddo
        enddo
        
        close(10)
    end subroutine
    
    !================================================================================================
    ! FUNCAO PARA LER UMA MATRIZ EM UM ARQUIVO .TXT NO ESQUEMA DE COORDENADAS
    !================================================================================================
    subroutine ReadMatCoord(A)
        implicit none
        
        type(RowPacked) :: A
        integer :: n, Density, i, indi, indj, tau
        integer, allocatable :: aux(:)
        real :: val
        
        Open(Unit = 10, File = "bcsstk01.txt", ACCESS = "SEQUENTIAL")
        read(10, *) n, tau
        
        Density = int(real(n*n)*MatrixDensity)+1
        allocate(aux(n),A%Len_Row(n),A%Row_Start(n),A%Col_Index(Density),A%Value(Density))
        
        A%Len_Row = 0; A%n = n
        
        do i = 1, tau
            read(10,*) indi, indj, val
            A%Len_Row(indi) = A%Len_Row(indi) + 1
        enddo
        
        close(10); Open(Unit = 10, File = "bcsstk01.txt", ACCESS = "SEQUENTIAL")
        read(10, *) n, tau
        
        A%Row_Start(1) = 1; aux(1) = 1
 
        do i = 2, n
            A%Row_Start(i) = A%Row_Start(i-1) + A%Len_Row(i-1)
            aux(i) = A%Row_Start(i)
        enddo

        do i = 1, tau
            read(10,*) indi, indj, val
            A%Col_Index(aux(indi)) = indj
            A%Value(aux(indi)) = Val
            aux(indi) = aux(indi) + 1
        enddo
        
        close(10)
    end subroutine ReadMatCoord
    
    !================================================================================================
    ! FUNCAO PARA LER UMA MATRIZ EM UM ARQUIVO .TXT NO ESQUEMA DE COORDENADAS
    !================================================================================================
    subroutine ReadMatCoord2(A)
        implicit none
        
        type(CoordPacked) :: A
        integer :: n, Density, i, indi, indj, tau
        integer, allocatable :: aux(:)
        real :: val
        
        Open(Unit = 10, File = "bcsstk01.txt", ACCESS = "SEQUENTIAL")
        read(10, *) n, tau
        
        allocate(A%Row_Index(tau),A%Col_Index(tau),A%Value(tau))
        
        A%nnz = n
        
        do i = 1, tau
            read(10,*) indi, indj, val
            A%Value = val
            A%Row_Index = indi
            A%Col_Index = indj
        enddo
        
        close(10)
    end subroutine ReadMatCoord2
    
end module
