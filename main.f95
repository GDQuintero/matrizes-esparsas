program main
    use gustavo
    use daniel
    
    real, dimension(5,5) :: A, B
    real, dimension(9,9) :: F, C
    type(RowPacked) :: D, E
    type(GustavsonPacked) :: G
    type(PivotMD) :: Pivo
    real, dimension(9) :: w
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
    
    F(1,:) = (/1., 2., 3., 0., 0., 0., 0., 0., 0./)
    F(2,:) = (/5., 6., 7., 8., 0., 0., 0., 0., 0./)
    F(3,:) = (/9., 1., 2., 3., 0., 0., 0., 0., 0./)
    F(4,:) = (/4., 5., 6., 7., 8., 0., 0., 0., 0./)
    F(5,:) = (/0., 0., 0., 9., 11., 0., 0., 0., 0./)
    F(6,:) = (/0., 0., 0., 0., 3., 4., 5., 6., 7./)
    F(7,:) = (/0., 0., 0., 0., 0., 8., 9., 1., 2./)
    F(8,:) = (/0., 0., 0., 0., 0., 3., 4., 5., 6./)
    F(9,:) = (/0., 0., 0., 0., 0., 7., 8., 9., 1./)
    
    
!     E = GatherRow(F)
    G = GatherGustavson(A)
    
!     C = Unpaking(E)
    
!     do i = 1, 9
!         print*, C(i,:)
!     enddo
!     print* 
!     print*, E%Len_Row
!     call GaussElimination(E,p)
    

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
                
                
            enddo           
        else
            call system("clear")
            print*, "Ainda nao foi implementado!!!"
            return
        endif
        

    end subroutine GaussElimination
    
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
