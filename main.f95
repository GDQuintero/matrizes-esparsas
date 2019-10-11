program main
    use gustavo
    use daniel
    
    real, dimension(5,5) :: A, B, H
    real, dimension(9,9) :: F, C
    type(RowPacked) :: D, E
    type(GustavsonPacked) :: G
    type(PivotMD) :: Pivo
    real, dimension(6) :: w
    integer, dimension(9) :: P
    integer :: NonZero = 0
    w = 0.; P = (/1, 2, 3, 4, 5, 6, 7, 8, 9/)
    
    A(1,:) = (/1.d0, 0.d0, 0.d0, -1.d0, 0.d0/)
    A(2,:) = (/2.d0, 0.d0, -2.d0, 0.d0, 3.d0/)
    A(3,:) = (/0.d0, -3.d0, 0.d0, 0.d0, 0.d0/)
    A(4,:) = (/0.d0, 4.d0, 0.d0, -4.d0, 0.d0/)
    A(5,:) = (/5.d0, 0.d0, -5.d0, 0.d0, 6.d0/)
    
    B(1,:) = (/1.d0, 0.d0, 0.d0, -1.d0, 5.d0/)
    B(2,:) = (/0.d0, 2.d0, -2.d0, 0.d0, 3.d0/)
    B(3,:) = (/0.d0, -2.d0, 1.d0, 0.d0, 0.d0/)
    B(4,:) = (/-1.d0, 0.d0, 0.d0, -4.d0, 0.d0/)
    B(5,:) = (/5.d0, 3.d0, 0.d0, 0.d0, 6.d0/)
    
    F(1,:) = (/1., 2., 3., 4., 0., 0., 0., 0., 0./)
    F(2,:) = (/5., 6., 7., 8., 0., 0., 0., 0., 0./)
    F(3,:) = (/9., 1., 2., 3., 0., 0., 0., 0., 0./)
    F(4,:) = (/4., 5., 6., 7., 8., 0., 0., 0., 0./)
    F(5,:) = (/0., 0., 0., 9., 1., 2., 0., 0., 0./)
    F(6,:) = (/0., 0., 0., 0., 3., 4., 5., 6., 7./)
    F(7,:) = (/0., 0., 0., 0., 0., 8., 9., 1., 2./)
    F(8,:) = (/0., 0., 0., 0., 0., 3., 4., 5., 6./)
    F(9,:) = (/0., 0., 0., 0., 0., 7., 8., 9., 1./)
    
    
    E = GatherRow(F)
    
    call GaussElimination(E,p)
    

    contains

    !================================================================================================
    ! UM PASSO DA ELIMINACAO DE GAUSS USANDO PIVOTAMENTO LOCAL
    !================================================================================================
    subroutine GaussElimination(A,P)
        implicit none
        
        type(RowPacked) :: A
        type(PivotMD) :: Pivo
        integer :: P(:), Criterio, i, j, k, n, ind = 0
        integer, allocatable :: tmp(:)
        real :: Mult=1.d0
        
        call system("clear")
        n = size(A%Len_Row); allocate(tmp(size(p)))
        print*, "Escolha uma estrategia de Pivotamento Local (digite apenas o numero): "
        print*, "1: Grau minimo"
        print*, "2: Minimum Fill-in"
        read*, Criterio
        
        if (Criterio .eq. 1) then
            call system("clear")
            
            do i = 1, n-1
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
                    
                enddo
            enddo           
        else
            call system("clear")
            print*, "Ainda nao foi implementado!!!"
            return
        endif
        

    end subroutine GaussElimination
    
end program
