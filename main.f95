program main
    use gustavo
    use daniel
    
    real, dimension(5,5) :: A, B, C
    type(RowPacked) :: D, E
    type(Pivot) :: Pivo
    real, dimension(5) :: w
    w = 0.d0
    
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
    
    E = GatherRow(B)
!     D = OneStepGaussElimination(E)
    call row_perm_rowpacked(E,1,3)
    call col_perm_rowpacked(E,1,3)
    call RowSumColPacked(E,2,1,1.,w)

    C = Unpaking(E)
    do i = 1, 5
        print*, C(i,:)
    enddo
    
    contains

    !================================================================================================
    ! UM PASSO DA ELIMINACAO DE GAUSS USANDO PIVOTAMENTO LOCAL
    !================================================================================================
    function OneStepGaussElimination(A)
        implicit none
        
        type(RowPacked) :: OneStepGaussElimination, A
        type(Pivot) :: Pivo
        integer :: Criterio, i, j, n
        real :: Mult=1.d0, ValPivo = 0
        
        call system("clear")
        n = size(A%Len_Row)
        print*, "Escolha um Critério de Pivotamento Local (digite apenas o numero): "
        print*, "1: Criterio de Markowitz"
        print*, "2: Criterio de Grau Minimo"
        read*, Criterio
        
        if (Criterio .eq. 1) then
            print*, "Ate que enfim"
            return
        elseif (Criterio .eq. 2) then
            call system("clear")
            Pivo = MinDeg(A)
            do i = 1, A%Len_Row(Pivo%Row)
                if (A%Col_Index(A%Row_Start(Pivo%Row)+i-1) .eq. Pivo%Col) then
                    ValPivo = A%Value(A%Row_Start(Pivo%Row)+i-1)
                    exit
                endif
            enddo
                        
            call row_perm_rowpacked(A,1,Pivo%Row)
            call col_perm_rowpacked(A,1,Pivo%Col)
            
            do i = 2, n
                do j = 1, A%Len_Row(i)
                    if (A%Col_Index(A%Row_Start(i)+j-1) .eq. ValPivo) then
                        Mult = ValPivo / A%Value(A%Row_Start(i)+j-1)
                        call RowSumColPacked(A,i,1,-Mult,w)
                    endif                        
                enddo
            enddo
            return
            
        else
            print*, "Erro: Digitou uma opcao invalida"
            return
        endif
        
!         call col_permutation(A,1,Pivo%Col)
    end function OneStepGaussElimination
    
end program
