!  XDR Fortran Interface
!  2014 (c) James W. Barnett <jbarnet4@tulane.edu> + SoX
!  https://github.com/wesbarnett/

module xtc

    use, intrinsic :: iso_c_binding, only: C_PTR, C_CHAR, C_FLOAT, C_INT

    implicit none
    private
    public xdrfile_open, xdrfile_close, read_xtc_natoms, read_xtc, write_xtc, read_xtc_n_frames

    ! the data type located in libxdrfile
    type, public, bind(C) :: xdrfile
      type(C_PTR) :: fp, xdr
      character(kind=C_CHAR) :: mode
      integer(C_INT) :: buf1, buf1size, buf2, buf2size
    end type xdrfile

    ! interface with libxdrfile
    interface 

      integer(C_INT) function read_xtc_natoms(filename,NATOMS) bind(C, name='read_xtc_natoms')
        import
        character(kind=C_CHAR), intent(in) :: filename
        integer(C_INT), intent(out) :: NATOMS
      end function

      type(C_PTR) function xdrfile_open(filename,mode) bind(C, name='xdrfile_open')
        import
        character(kind=C_CHAR), intent(in) :: filename(*), mode(*)
      end function

      integer(C_INT) function read_xtc(xd,NATOMS,STEP,time,box,x,prec) bind(C, name='read_xtc')
        import
        type(xdrfile), intent(in) :: xd
        integer(C_INT), intent(out) :: NATOMS, STEP
        real(C_FLOAT), intent(out) :: time, prec, box(*), x(*)
      end function

      integer(C_INT) function xdrfile_close(xd) bind(C,name='xdrfile_close')
        import
        type(xdrfile), intent(in) :: xd
      end function

      integer(C_INT) function write_xtc(xd,NATOMS,STEP,time,box,x,prec) bind(C)
        import
        type(xdrfile), intent(in) :: xd
        integer(C_INT), intent(in), value :: NATOMS, STEP
        real(C_FLOAT), intent(in), value :: time, prec
        real(C_FLOAT), intent(in) :: box(*), x(*)
      end function

      integer(C_INT) function read_xtc_n_frames(filename, N_FRAMES, EST_NFRAMES, OFFSETS) bind(C,name="read_xtc_n_frames")
        import
        character(kind=C_CHAR), intent(in) :: filename(*)
        integer(C_INT), intent(out) :: N_FRAMES, EST_NFRAMES
        type(C_PTR) :: OFFSETS
      end function

    end interface

end module xtc