program main
    use gustavo
!     use daniel
    
    real, dimension(5,5) :: A, B
    real, dimension(5) :: x, y, w
    type(PackedVector) :: xx, yy, soma
    type(RowPacked) :: C, E
    type(ColPacked) :: D
    
    w = 0.d0
    
    A(1,:) = (/1.d0, 0.d0, 0.d0, -1.d0, 0.d0/)
    A(2,:) = (/2.d0, 0.d0, -2.d0, 0.d0, 3.d0/)
    A(3,:) = (/0.d0, -3.d0, 0.d0, 0.d0, 0.d0/)
    A(4,:) = (/0.d0, 4.d0, 0.d0, -4.d0, 0.d0/)
    A(5,:) = (/5.d0, 0.d0, -5.d0, 0.d0, 6.d0/)
    
    
    C = GatherRow(A)
    
end program
