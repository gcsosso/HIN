module MOD_vector3
    
    type :: vector3
        integer :: rings(3)
    end type vector3
    
    type :: vector4
        integer :: rings(4)
    end type vector4
    
    type :: vector_alloc
        integer, allocatable :: rings(:)
    end type vector_alloc
    
    type :: vector
        integer, dimension(:,:), allocatable :: mrings
    end type vector
    type :: ragged_array
        type(vector), dimension(:), allocatable :: stat_wr_size
    end type ragged_array
    
    type :: vector2
        integer :: atom_match(2)
    end type vector2
    type :: cnx
        type(vector2) :: matches(5)
    end type cnx
    type :: cnx_graph
        type(cnx), dimension(:), allocatable :: ring_cnx
    end type cnx_graph
    
end module