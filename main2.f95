program main2
	use gustavo
	implicit none


	contains

	!================================================================================================
    ! LEER UNA MATRIZ DENSA Y ALMACENARLA EN EL FORMATO GUSTAVSON
    !================================================================================================
    subroutine ReadMatGus(A)
    	implicit none

    	type(GustavsonPacked) :: A
        integer :: n, Density, i, j, k, NonZero
        real, allocatable :: Numbers(:)
        
        Open(Unit = 10, File = "dense1.txt", ACCESS = "SEQUENTIAL")
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
    end subroutine
end program