
#include "math_def.h"
Module mesh
 Implicit None

 Type t_vertex
  Integer :: id !this vertex' position inside the 'vertices' array
  Integer :: submesh_id !the submesh this vertex belongs to
  Real (Kind=DEF_DBL_PREC) :: x, y, z
  Complex (Kind=DEF_DBL_PREC), Allocatable :: vec(:) !wave function at this point on the fermi surface
  Type (t_vertex), Allocatable :: nbs(:) !list of ids of nbs
 End Type

 Type t_mesh
  Integer :: alloc = 0
  Type (t_vertex), Allocatable :: vertices(:)
 End Type

Contains
 
 Subroutine alloc_mesh(mesh,n)
  Implicit None
  Type (t_mesh) :: mesh
  Integer :: n
  !local
  
  If (mesh%alloc == 1) Return ! cannot allocate previously allocated mesh 

  Allocate(mesh%vertices(n))
  mesh%alloc = 1
 
  Return
 End Subroutine
  
 Subroutine free_mesh(mesh)
  Implicit None
  Type (t_mesh) :: mesh
  !local
  Integer :: i

  if (mesh%alloc == 0) Return ! cannot free non-allocated mesh

  Do i=1,size(mesh%vertices)
   Deallocate(mesh%vertices(i)%nbs)
   Deallocate(mesh%vertices(i)%vec)
  Enddo
 
  Deallocate(mesh%vertices)
  mesh%alloc = 0

  Return
 End Subroutine
End Module
