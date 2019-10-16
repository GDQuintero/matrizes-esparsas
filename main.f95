program main
    use gustavo
    use daniel
    
    type(RowPacked) :: A
    type(PivotMD) :: Pivo
    real, allocatable :: w(:), B(:,:), per(:,:), C(:,:), D(:,:)
    character(len=3), allocatable :: E(:,:)
    integer, allocatable :: P(:)
    
    call ReadMatrix(A)
    allocate(w(A%n),P(A%n),B(A%n,A%n),per(A%n,A%n),C(A%n,A%n),D(A%n,A%n))
    
    do i = 1, A%n
        P(i) = i
    enddo
    
    call GaussElimination(A,P)
    call Export(A,A%n)
    B = Unpaking(A)
    per = PerMat(P)
    C = ProdMat(transpose(per),B); D = ProdMat(C,per)
    
!     
!     call PrintMat(D,A%n)
!     print*
!     E = Pattern(D,A%n)
    
!     print*, p
    contains
    
    !================================================================================================
    ! CRITERIO DE GRADO MINIMO
    !================================================================================================
    function MinDeg(A,P,ind)
        implicit none
        
        integer :: i = 0, r_min = 0
        type(RowPacked) :: A
        type(PivotMD) :: MinDeg
        integer :: P(:), ind
        
        !Escojemos la primera fila no ignorada que aparece en el vector P
        r_min = A%Len_Row(P(ind+1))
        
        !Verificamos si la primera fila considerada es singleton
        if (r_min .eq. 1) then
            MinDeg%Row = P(ind+1)
            MinDeg%Col = MinDeg%Row
            MinDeg%Value = A%Value(A%Row_Start(P(ind+1)))
            MinDeg%rmin = r_min
            return
        endif
        
        MinDeg%Row = P(ind+1)
        MinDeg%Col = MinDeg%Row
        MinDeg%rmin = r_min
        
        !Buscamos min r_i, y seleccionmos el elemento diagonal como pivote
        do i = ind+2, size(A%Len_Row)
            !Fila singleton
            if (A%Len_Row(P(i)) .eq. 1) then
                r_min = A%Len_Row(P(i))
                MinDeg%Row = P(i)
                MinDeg%Col = MinDeg%Row
                MinDeg%Value = A%Value(A%Row_Start(P(i)))
                MinDeg%rmin = r_min
                return
            elseif (A%Len_Row(P(i)) .lt. r_min) then
                r_min = A%Len_Row(P(i))
                MinDeg%Row = P(i)
                MinDeg%Col = MinDeg%Row
                MinDeg%rmin = r_min
            endif
        enddo
        
        !Buscamos el valor del pivote
        do i = 1, r_min
            if (A%Col_Index(A%Row_Start(MinDeg%Row)+i-1) .eq. MinDeg%Row) then
                MinDeg%Value = A%Value(A%Row_Start(MinDeg%Row)+i-1)
            endif
        enddo
        
    end function MinDeg
    
    !================================================================================================
    ! ELIMINACAO DE GAUSS USANDO PIVOTAMENTO LOCAL
    !================================================================================================
    subroutine GaussElimination(A,P)
        implicit none
        
        type(RowPacked) :: A
        type(PivotMD) :: Pivo
        integer :: P(:), Criterio, i=0, j=0, k=0, n=0, ind = 0
        integer, allocatable :: tmp(:)
        real :: Mult=1.d0
        
!         call system("clear")
        n = A%n; allocate(tmp(size(p)))
        
!         call system("clear")
        
        do i = 1, n - 1
            Pivo = MinDeg(A,P,ind)!Calculamos el pivote
            ind = ind + 1!Indice para ignorar la fila de los pivotes elegidos
            tmp = P
            P(ind) = Pivo%Row
            
            do j = ind + 1, n
                if (P(j) .eq. Pivo%Row) then
                    P(j) = tmp(ind)
                    exit!Hasta aqui solo para guardar la permutacion
                endif
            enddo

            do j = ind + 1, n
                do k = A%Row_Start(P(j)), A%Row_Start(P(j)) + A%Len_Row(P(j))-1
                    if (A%Col_Index(k) .eq. Pivo%Row) then
                        Mult = -1.d0*A%Value(k) / Pivo%Value
                        call RRowSumRowPacked(A,P(j),Pivo%Row,Mult,w)
                    endif
                enddo
            enddo
        enddo           
        
        
    end subroutine GaussElimination
    
    !================================================================================================
    ! FUNCAO PARA LER UMA MATRIZ EM UM ARQUIVO .TXT
    !================================================================================================
    subroutine ReadMatrix(A)
        implicit none
        
        type(RowPacked) :: A
        integer :: n, Density, i, j, k, NonZero
        real, allocatable :: Numbers(:)
        
        Open(Unit = 10, File = "matriz2.txt", ACCESS = "SEQUENTIAL")
        read(10, *) n
        
        A%n = n; A%m = n
        Density = int(real(n**2)*MatrixDensity)+1; NonZero = 0; k = 0
        allocate(Numbers(n),A%Len_Row(n),A%Row_Start(n),A%Col_Index(Density),A%Value(Density))
        
        do i = 1, n
            read(10, *) Numbers
            A%Row_Start(i) = k + 1
            do j = 1, n
                if (abs(Numbers(j)) .gt. 10D-4) then
                    k = k + 1
                    A%Col_Index(k) = j 
                    A%Value(k) = Numbers(j)
                    NonZero = NonZero + 1
                endif
            enddo
            A%Len_Row(i) = NonZero
            NonZero = 0
        enddo
        close(10)
    end subroutine ReadMatrix
    
    !================================================================================================
    ! FUNCAO PARA LER UMA MATRIZ EM UM ARQUIVO .TXT NO ESQUEMA DE COORDENADAS
    !================================================================================================
    subroutine ReadMatCoord(A)
        implicit none
        
        type(RowPacked) :: A
        integer :: n, Density, i, indi, indj, tau
        integer, allocatable :: aux(:)
        real :: val
        
        Open(Unit = 10, File = "bcsstm07.txt", ACCESS = "SEQUENTIAL")
        read(10, *) n, tau
        
        Density = int(real(n*n)*MatrixDensity)+1
        allocate(aux(n),A%Len_Row(n),A%Row_Start(n),A%Col_Index(Density),A%Value(Density))
        
        A%Len_Row = 0; A%n = n
        
        do i = 1, tau
            read(10,*) indi, indj, val
            A%Len_Row(indi) = A%Len_Row(indi) + 1
        enddo
        close(10); Open(Unit = 10, File = "bcsstm07.txt", ACCESS = "SEQUENTIAL")
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
end program
