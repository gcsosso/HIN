module MOD_vector3
    
    type :: vector3
        integer :: rings(3)
    end type vector3
    
    type :: vector
        integer, dimension(:,:), allocatable :: mrings
    end type vector
    type :: ragged_array
        type(vector), dimension(:), allocatable :: stat_wr_size
    end type ragged_array
    
end module