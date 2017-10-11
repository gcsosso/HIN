module dfs
  integer,allocatable,dimension(:,:) :: graph_solid_connect, lwho
  integer,allocatable,dimension(:) :: neigh, followgraph, predecessor
  integer ,allocatable,dimension(:) :: volume, cr_list
  integer :: black, count_cls
end module dfs
