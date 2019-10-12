program main
    use gustavo
    use daniel
    
    type(RowPacked) :: A
    type(PivotMD) :: Pivo
    real, allocatable :: w(:), B(:,:)
    integer, allocatable :: P(:)
    integer :: NonZero = 0
    
    call ReadMatrix(A)
    allocate(w(A%n),P(A%n),B(A%n,A%n))
    
    do i = 1, A%n
        P(i) = i
    enddo
    
    call GaussElimination(A,P)
    
    B = Unpaking(A)
    
    do i = 1, A%n
        print*, B(i,:)
    enddo
    
    print*
!     print*, "As permutacoes sao: ", p
    contains

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
        
        call system("clear")
        n = A%n; allocate(tmp(size(p)))
        print*, "Escolha uma estrategia de Pivotamento Local (digite apenas o numero): "
        print*, "1: Grau minimo"
        print*, "2: Minimum Fill-in"
        read*, Criterio
        
        if (Criterio .eq. 1) then
            call system("clear")
            
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
        else
            call system("clear")
            print*, "Ainda nao foi implementado!!!"
            return
        endif
        
    end subroutine GaussElimination
    
end program
