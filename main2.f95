program main2
    use gustavo
    implicit none

    type(GustavsonPacked) :: A

    call ReadMatGus(A)
    print*, A%Len_Row
	contains

    !================================================================================================
    ! LEER UNA MATRIZ DENSA Y ALMACENARLA EN EL FORMATO GUSTAVSON - ORDEM N*TAU
    !================================================================================================
    subroutine ReadMatGus(A)
    	implicit none

        type(GustavsonPacked) :: A
        integer :: n, nz, Density, i, j, rk, ck, Col_nz, Row_nz
        real(kind=8), allocatable :: Numbers(:), aux(:,:)
        real(kind=8) :: tol = 10d-4
        
        Open(Unit = 10, File = "dense2.txt", ACCESS = "SEQUENTIAL")
        read(10, *) n, nz

        A%n = n; A%nz = nz; Col_nz = 0; Row_nz = 0
        rk = 0; ck = 0; Density = int(real(n**2)*MatrixDensity)+1 
        allocate(Numbers(n),A%Len_Row(n),A%Len_Col(n),A%Row_Start(n),A%Col_Start(n))
        allocate(A%Value(Density),A%Row_Index(Density),A%Col_Index(Density),aux(n,n))
        A%Value = 0.d0; A%Len_Col = 0; A%Len_Row = 0; A%Col_Start = 0; A%Row_Start = 0
        A%Row_Index = 0; A%Col_Index = 0; aux = 0.d0

        do i = 1, n
            read(10, *) Numbers
            A%Row_Start(i) = rk + 1
          
            do j = 1, n
                if (abs(Numbers(j)) .gt. tol) then
                    aux(i,j) = Numbers(j)
                    rk = rk + 1
                    A%Col_Index(rk) = j 
                    Row_nz = Row_nz + 1
                endif
            enddo
            A%Len_Row(i) = Row_nz
            Row_nz = 0
        enddo

        do j = 1, n
            A%Col_Start(j) = ck + 1
            
        	do i = 1, n
        		if (abs(aux(i,j)) .gt. tol) then
                    ck = ck + 1
                    A%Row_Index(ck) = i
                    A%Value(ck) = aux(i,j)
                    Col_nz = Col_nz + 1
        		endif
        	enddo
            A%Len_Col(j) = Col_nz
            Col_nz = 0
        enddo
        close(10)

    end subroutine
end program