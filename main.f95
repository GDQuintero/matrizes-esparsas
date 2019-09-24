program main
    use gustavo
    use daniel
    
    real, dimension(5,5) :: A, B
    real, dimension(9,9) :: F, C
    type(RowPacked) :: D, E
    type(ColPacked) :: G
    type(Pivot) :: Pivo
    real, dimension(5) :: w
    integer :: NonZero = 0
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
    
    F(1,:) = (/1., 2., 3., 4., 0., 0., 0., 0., 0./)
    F(2,:) = (/5., 6., 7., 8., 0., 0., 0., 0., 0./)
    F(3,:) = (/9., 1., 2., 3., 0., 0., 0., 0., 0./)
    F(4,:) = (/4., 5., 6., 7., 8., 0., 0., 0., 0./)
    F(5,:) = (/0., 0., 0., 9., 4., 2., 0., 0., 0./)
    F(6,:) = (/0., 0., 0., 0., 3., 4., 5., 6., 7./)
    F(7,:) = (/0., 0., 0., 0., 0., 8., 9., 1., 2./)
    F(8,:) = (/0., 0., 0., 0., 0., 3., 4., 5., 6./)
    F(9,:) = (/0., 0., 0., 0., 0., 7., 8., 9., 1./)
    
    
    E = GatherRow(B)
    call OneStepGaussElimination(E)
    C = Unpaking(E)
    
    do i = 1, 5
        print*, C(i,:)        
    enddo
    
!     do i = 1, 9
!         do j = 1, 9
!             if (C(i,j) .ne. 0) then
!                 NonZero = NonZero + 1
!             endif
!         enddo
!     enddo
!     print*
!     print*, NonZero

!     G = GatherCol(A)
!     Pivo = MinFillin(G)
    contains

    !================================================================================================
    ! UM PASSO DA ELIMINACAO DE GAUSS USANDO PIVOTAMENTO LOCAL
    !================================================================================================
    subroutine OneStepGaussElimination(A)
        implicit none
        
        type(RowPacked) :: A
        type(ColPacked) :: B
        type(Pivot) :: Pivo
        integer :: Criterio, i, j, n
        real :: Mult=1.d0
        
        call system("clear")
        n = size(A%Len_Row)
        print*, "Escolha uma estrategia de Pivotamento Local (digite apenas o numero): "
        print*, "1: Criterio de Markowitz"
        print*, "2: Grau minimo"
        print*, "3: Min. row in min. col"
        read*, Criterio
        
        if (Criterio .eq. 1) then
            call system("clear")
            B = PackRowCol(A)
            Pivo = Markowitz(B)
            
        elseif (Criterio .eq. 2) then
            call system("clear")
            Pivo = MinDeg(A)
        
        elseif (Criterio .eq. 3) then
            call system("clear")
            B = PackRowCol(A)
            Pivo = MinRowInMinCol(B)
        
        else
            print*, "Erro: Digitou uma opcao invalida"
            return
        endif
        print*, pivo%row , pivo%col, Pivo%Value
        print*
!         do i = 1, A%Len_Row(Pivo%Row)
!             if (A%Col_Index(A%Row_Start(Pivo%Row)+i-1) .eq. Pivo%Col) then
!                 ValPivo = A%Value(A%Row_Start(Pivo%Row)+i-1)
!                 exit
!             endif
!         enddo
                    
        call row_perm_rowpacked(A,1,Pivo%Row)
        call col_perm_rowpacked(A,1,Pivo%Col)
        
        do i = 2, n
            do j = 1, A%Len_Row(i)
                if (A%Col_Index(A%Row_Start(i)+j-1) .eq. Pivo%Value) then
                    Mult = Pivo%Value / A%Value(A%Row_Start(i)+j-1)
                    call RowSumColPacked(A,i,1,-Mult,w)
                endif                        
            enddo
        enddo
        
    end subroutine OneStepGaussElimination
    
    !================================================================================================
    !  MIN FILL-IN (Forma ColPacked)
    !================================================================================================
    function MinFillin(A)
        implicit none
        
        integer :: i, j, k=0, n, l=0
        type(ColPacked) :: A
        type(RowPacked) :: aux
        type(Pivot) :: MinFillin
        type(EntryPacked) :: Pivos
        real :: Mult
        
        n = size(A%Len_Col); aux = PackColRow(A)
        
        !Alocamos o type contendo os possiveis pivos
        allocate(Pivos%Row_Index(n*n),Pivos%Col_Index(n*n),Pivos%Value(n*n))
        
        !Procuramos os pivos candidatos
        do i = 1, n
            do j = A%Col_Start(i), A%Col_Start(i) + A%Len_Col(i) - 1
                if (abs(A%Value(j)) .ge. u*maxval(abs(A%Value(A%Col_Start(i):A%Col_Start(i) + A%Len_Col(i) - 1)))) then
                    l = l + 1
                    Pivos%Row_Index(l) = A%Row_Index(j)
                    Pivos%Col_Index(l) = i
                    Pivos%Value(l) = A%Value(j)
                endif
            enddo
        enddo
!         print*, Pivos%Row_Index(1:l)
!         print*, Pivos%Col_Index(1:l)
!         print*, Pivos%Value(1:l)
        
        do k = 1, 1
            call row_perm_rowpacked(aux,1,Pivos%Row_Index(k))
            call col_perm_rowpacked(aux,1,Pivos%Col_Index(k))
            
            do i = 2, n
                print*, aux%Len_Row(i)
                do j = 1, aux%Len_Row(i)
                    if (aux%Col_Index(aux%Row_Start(i)+j-1) .eq. Pivos%Value(k)) then
                        Mult = Pivos%Value(k) / aux%Value(aux%Row_Start(i)+j-1)
!                         call RowSumColPacked(aux,i,1,-Mult,w)
                    endif                        
                enddo
            enddo
!             print*, aux%Len_Row
        enddo
    end function MinFillin
    
end program
