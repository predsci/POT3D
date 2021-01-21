c
c-----------------------------------------------------------------------
c
c ****** Source to build the parsing library.
c ****** These routines are used by Zoran Mikic's tools.
c
c **********************************************************************
c
c Copyright 2018 Predictive Science Inc.
c
c Licensed under the Apache License, Version 2.0 (the "License");
c you may not use this file except in compliance with the License.
c You may obtain a copy of the License at
c
c    http://www.apache.org/licenses/LICENSE-2.0
c
c Unless required by applicable law or agreed to in writing, software
c distributed under the License is distributed on an "AS IS" BASIS,
c WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
c implied.
c See the License for the specific language governing permissions and
c limitations under the License.
c
c **********************************************************************
c
c
c-----------------------------------------------------------------------
c
c        07/29/2003, ZM, Version 1.00:
c
c         - Original version of the parsing library.
c           This library was put together to facilitate the
c           development of ZM's tools.
c           It includes routines to parse the command line.
c           The code was cleaned up to use standard FORTRAN90.
c
c        01/17/2005, ZM, Version 1.01:
c
c         - Added the function NARGS_SPECIFIED to return the
c           number of arguments specified on the command line.
c
c        10/28/2005, ZM, Version 1.02:
c
c         - Added the functions LOAD_LIST_OF_REALS and
c           LOAD_LIST_OF_INTS that can be used to parse
c           arguments that contain lists of real or integer values.
c         - Changed the length of temporary arguments to equal
c           512 characters to accomodate long file names.
c
c        10/31/2006, ZM, Version 1.03:
c
c         - Removed the EXTERNAL declarations for GETARG and IARGC
c           since these are now intrinsic routines in the
c           Intel 9.1 compiler.
c
c        03/10/2008, ZM, Version 1.04:
c
c         - Added the LCASE and UCASE functions to convert strings
c           to lowercase and uppercase.
c
c        04/27/2018, ZM, Version 1.05:
c
c         - Added the ability to specify group-sets of keywords.
c           This extends the flexibility of the syntax allowed
c           for keywords.  For example, it allows a syntax in which
c           only one of a set of keywords is allowed to be specified
c           at a time.
c
c        01/08/2019, RC, Version 1.06:
c
c         - Added "append" option to ffopen.
c
c-----------------------------------------------------------------------
c
c#######################################################################
      module parselib_ident
c
      character(*), parameter :: cname='PARSELIB'
      character(*), parameter :: cvers='1.06'
      character(*), parameter :: cdate='01/08/2019'
c
      end module
c#######################################################################
      module parse_args
c
      use string_def
c
      implicit none
c
c-----------------------------------------------------------------------
c ****** Argument descriptor and storage for command-line arguments.
c-----------------------------------------------------------------------
c
c ****** Structure to hold a group-set.
c
      type :: group_set
        integer :: group_type
        integer :: n_members
        integer, dimension(:), allocatable :: members
        logical :: required
        logical :: only_one
      end type
c
c ****** Structure to hold an argument.
c
      type :: arg_descriptor
        integer :: group
        logical :: set
        logical :: required
        type(string) :: keyword
        type(string) :: name
        type(string) :: value
      end type
c
c ****** Maximum number of arguments.
c
      integer, parameter :: mxarg=100
c
c ****** Number of arguments defined.
c
      integer :: nargs
c
c ****** Argument descriptor.
c
      type(arg_descriptor), dimension(mxarg) :: args
c
c ****** Number of arguments specified.
c
      integer :: nargs_spec
c
c ****** Maximum number of group-sets.
c
      integer, parameter :: mxgroup_sets=10
c
c ****** Number of group-sets defined.
c
      integer :: ngroup_sets=0
c
c ****** Group-set descriptor.
c
      type(group_set), dimension(mxgroup_sets) :: group_sets
c
      end module
c#######################################################################
      subroutine ffopen (iun,fname,mode,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Open file FNAME and link it to unit IUN.
c
c ****** If there is an error, this routine returns IERR.ne.0.
c
c-----------------------------------------------------------------------
c
c ****** When MODE='r', the file must exist.
c ****** When MODE='w', the file is created.
c ****** When MODE='rw', the file must exist, but can be overwritten.
c ****** When MODE='a', the file is created if it does not exist,
c ******                otherwise, it is appended.
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: iun
      character(*) :: fname
      character(*) :: mode
      integer :: ierr
      logical :: ex
c
c-----------------------------------------------------------------------
c
      ierr=0
c
      if (mode.eq.'r') then
        open (iun,file=fname,form="FORMATTED",status='old',err=900)
      else if (mode.eq.'rw') then
        open (iun,file=fname,form="FORMATTED",status='replace',err=900)
      else if (mode.eq.'w') then
        open (iun,file=fname,form="FORMATTED",status='new',err=900)
      elseif (mode.eq.'a') then
        inquire(file=fname, exist=ex)
        if (ex) then
          open (iun,file=fname,form="FORMATTED",
     &          position='append',err=900)
        else
          open (iun,file=fname,form="FORMATTED",status='new',err=900)
        end if
      else
        write (*,*)
        write (*,*) '### ERROR in FFOPEN:'
        write (*,*) '### Invalid MODE requested.'
        write (*,*) 'MODE = ',mode
        write (*,*) 'File name: ',trim(fname)
        ierr=2
        return
      end if
c
      return
c
  900 continue
c
      write (*,*)
      write (*,*) '### ERROR in FFOPEN:'
      write (*,*) '### Error while opening the requested file.'
      write (*,*) 'File name: ',trim(fname)
      write (*,*) 'MODE = ',mode
      ierr=1
c
      return
      end
c#######################################################################
      function lcase (s)
c
c-----------------------------------------------------------------------
c
c ****** Convert the string S into lowercase letters and return it as
c ****** the function result.
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(*), intent(in) :: s
      character(len(s)) :: lcase
c
c-----------------------------------------------------------------------
c
      integer :: i,ic
c
c-----------------------------------------------------------------------
c
      lcase=' '
c
      do i=1,len_trim(s)
        ic=iachar(s(i:i))
        if (ic.ge.65.and.ic.le.90) then
          ic=ic+32
        end if
        lcase(i:i)=achar(ic)
      end do
c
      return
      end
c#######################################################################
      function ucase (s)
c
c-----------------------------------------------------------------------
c
c ****** Convert the string S into uppercase letters and return it as
c ****** the function result.
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(*), intent(in) :: s
      character(len(s)) :: ucase
c
c-----------------------------------------------------------------------
c
      integer :: i,ic
c
c-----------------------------------------------------------------------
c
      ucase=' '
c
      do i=1,len_trim(s)
        ic=iachar(s(i:i))
        if (ic.ge.97.and.ic.le.122) then
          ic=ic-32
        end if
        ucase(i:i)=achar(ic)
      end do
c
      return
      end
c#######################################################################
      subroutine parse (errmsg,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Parse the command line.
c
c-----------------------------------------------------------------------
c
c ****** The syntax for the keyword/argument items can be defined
c ****** by using routine DEFARG.
c
c ****** On return, IERR=0 indicates that the command line was
c ****** parsed successfully.
c
c ****** IERR=1 indicates that no arguments were present.  This
c ****** is usually used to print the usage line.
c
c ****** IERR=2 indicates that a syntax error occured.
c
c ****** IERR=3 indicates that one or more required arguments
c ****** was not supplied.
c
c ****** When IERR=2 or IERR=3, an error message is put into
c ****** character string ERRMSG.
c
c-----------------------------------------------------------------------
c
      use syntax
      use parse_args
      use get_str_interface
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(*) :: errmsg
      integer :: ierr
c
c-----------------------------------------------------------------------
c
c ****** Command line arguments.
c
      integer :: iargc
c
c-----------------------------------------------------------------------
c
c ****** Carriage return.
c
      character(1), parameter :: cr=achar(10)
c
c-----------------------------------------------------------------------
c
      integer, external :: matchkw
      integer, external :: nextarg
c
c-----------------------------------------------------------------------
c
      character(512) :: arg
      integer :: na,ia,ia0,iarg,ls,i,j,n
c
c-----------------------------------------------------------------------
c
c ****** Initialization.
c
      ierr=0
      nargs_spec=0
      errmsg=' '
c
c ****** Get the number of command line arguments.
c
      na=iargc()
      if (na.eq.0) then
        ierr=1
        return
      end if
c
      ia=1
  200 continue
c
      ia0=ia
c
c ****** Process arguments with syntax: <kw> <arg>
c
      if (na-ia+1.ge.2) then
        call getarg (ia,arg)
        iarg=matchkw(GROUP_KA,trim(arg))
        if (iarg.gt.0) then
          if (.not.args(iarg)%set) then
            ia=ia+1
            call getarg (ia,arg)
            call delete_str (args(iarg)%value)
            call put_str (trim(arg),args(iarg)%value)
            args(iarg)%set=.true.
            ia=ia+1
            nargs_spec=nargs_spec+1
            go to 300
          end if
        end if
      end if
c
c ****** Process arguments with syntax: <kw> <arg> <arg>
c
      if (na-ia+1.ge.3) then
        call getarg (ia,arg)
        iarg=matchkw(GROUP_KAA,trim(arg))
        if (iarg.gt.0) then
          if (.not.args(iarg)%set) then
            ia=ia+1
            call getarg (ia,arg)
            ls=len_trim(arg)
            ls=ls+1
            arg(ls:ls)=' '
            ia=ia+1
            call getarg (ia,arg(ls+1:))
            call delete_str (args(iarg)%value)
            call put_str (trim(arg),args(iarg)%value)
            args(iarg)%set=.true.
            ia=ia+1
            nargs_spec=nargs_spec+1
            go to 300
          end if
        end if
      end if
c
c ****** Process arguments with syntax: <kw>
c
      if (na-ia+1.ge.1) then
        call getarg (ia,arg)
        iarg=matchkw(GROUP_K,trim(arg))
        if (iarg.gt.0) then
          if (.not.args(iarg)%set) then
            call delete_str (args(iarg)%value)
            call put_str (' ',args(iarg)%value)
            args(iarg)%set=.true.
            ia=ia+1
            nargs_spec=nargs_spec+1
            go to 300
          end if
        end if
      end if
c
c ****** Process arguments with syntax: <arg>
c
      if (na-ia+1.ge.1) then
        iarg=nextarg(GROUP_A)
        if (iarg.gt.0) then
          call getarg (ia,arg)
          call delete_str (args(iarg)%value)
          call put_str (trim(arg),args(iarg)%value)
          args(iarg)%set=.true.
          ia=ia+1
          nargs_spec=nargs_spec+1
          go to 300
        end if
      end if
c
  300 continue
c
c ****** Check that an argument was found.
c
      if (ia.eq.ia0) then
        ierr=2
        errmsg='### Syntax error.'
        return
      end if
c
c ****** Keep processing arguments until done.
c
      if (na-ia+1.gt.0) go to 200
c
c ****** Check that the required arguments were supplied.
c
      do i=1,nargs
        if (args(i)%required.and..not.args(i)%set) then
          ierr=3
          errmsg='### A required argument was not supplied.'
          return
        end if
      enddo
c
c ****** Check that the group-set keywords were specified
c ****** properly.
c
c ****** Check the "only one" and "required" group-set keywords.
c
      do i=1,ngroup_sets
        n=0
        do j=1,group_sets(i)%n_members
          iarg=group_sets(i)%members(j)
          if (args(iarg)%set) n=n+1
        enddo
        if (n.gt.1.and.group_sets(i)%only_one) then
          ierr=2
          errmsg='### Syntax error.'
          errmsg=trim(errmsg)//cr//cr
          errmsg=trim(errmsg)//' ### You can only specify'//
     &           ' one of the following keywords:'
          do j=1,group_sets(i)%n_members
            iarg=group_sets(i)%members(j)
            errmsg=trim(errmsg)//cr//' '//
     &             get_str(args(iarg)%keyword)
          enddo
          return
        else if (n.eq.0.and.group_sets(i)%required) then
          ierr=2
          errmsg='### Syntax error.'
          errmsg=trim(errmsg)//cr//cr
          errmsg=trim(errmsg)//' ### You must specify'//
     &           ' one of the following keywords:'
          do j=1,group_sets(i)%n_members
            iarg=group_sets(i)%members(j)
            errmsg=trim(errmsg)//cr//' '//
     &             get_str(args(iarg)%keyword)
          enddo
          return
        end if
      enddo
c
      return
      end
c#######################################################################
      subroutine get_usage_line (usage)
c
c-----------------------------------------------------------------------
c
c ****** Construct the usage line in paragraph USAGE.
c
c ****** Use routine PRINT_PAR to write the usage line.
c
c-----------------------------------------------------------------------
c
      use syntax
      use parse_args
      use paragraph_def
      use new_par_interface
      use add_line_interface
      use get_str_interface
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(paragraph), pointer :: usage
c
c-----------------------------------------------------------------------
c
c ****** Right margin for printing the usage line.
c
      integer, parameter :: rmargin=78
c
c-----------------------------------------------------------------------
c
      character(512) :: line
      integer :: iarg,n0
      type(paragraph), pointer :: current_par
c
      logical, dimension(:), allocatable :: done
      integer :: i,j,ix
      logical :: belongs_to_gs
      character(512) :: str
c
c-----------------------------------------------------------------------
c
c ****** Construct the usage line in USAGE.
c
      call new_par (usage)
      current_par=>usage
c
c ****** Start with the command name (as invoked).
c
      call getarg (0,line)
c
c ****** Allocate storage.
c
      allocate (done(nargs))
      done=.false.
c
      iarg=1
c
c ****** Add the arguments one by one.
c
      do while (iarg.le.nargs)
c
        if (done(iarg)) then
          iarg=iarg+1
          cycle
        end if
c
c ****** Add the syntax for the next argument to LINE.
c
        n0=len_trim(line)
c
c ****** Check if this argument belongs to a "only one" group-set.
c ****** If so, collect all members of the group-set.  Otherwise, just
c ****** process this argument.
c
        belongs_to_gs=.false.
        do i=1,ngroup_sets
          if (group_sets(i)%only_one) then
            do j=1,group_sets(i)%n_members
              if (group_sets(i)%members(j).eq.iarg) then
                ix=i
                belongs_to_gs=.true.
                exit
              end if
            enddo
          end if
          if (belongs_to_gs) exit
        enddo
c
        if (belongs_to_gs) then
c
c ****** This argument belongs to a "only one" group-set:
c ****** process all its members.
c
          str=''
          do j=1,group_sets(ix)%n_members
            i=group_sets(ix)%members(j)
            done(i)=.true.
            str=trim(str)//get_str(args(i)%keyword)
            if (j.lt.group_sets(ix)%n_members) then
              str=trim(str)//'|'
            end if
          enddo
          if (group_sets(ix)%required) then
            line=trim(line)//' '//trim(str)
          else
            line=trim(line)//' ['//trim(str)//']'
          end if
c
        else
c
c ****** This argument does not belong to a "only one" group-set:
c ****** process just this single argument.
c
          if (args(iarg)%required) then
            line=trim(line)//' '//get_str(args(iarg)%keyword)
          else
            line=trim(line)//' ['//get_str(args(iarg)%keyword)
          end if
          line=trim(line)//' '//get_str(args(iarg)%name)
          if (.not.args(iarg)%required) then
            line=trim(line)//']'
          end if
          done(iarg)=.true.
c
        end if
c
c ****** Check if the addition of the argument causes the line
c ****** to wrap; if it does, break the line prior to the
c ****** argument text.
c
        if (len_trim(line).gt.rmargin) then
          call add_line (line(1:n0),current_par)
          line=' '//line(n0+1:)
        end if
c
c ****** If the line is still too long, force a break at RMARGIN
c ****** until the line is shorter than RMARGIN.
c
        do while (len_trim(line).gt.rmargin)
          call add_line (line(1:rmargin),current_par)
          line='  '//line(rmargin+1:)
        enddo
c
c ****** Process the next argument.
c
        iarg=iarg+1
c
      enddo
c
c ****** Add the last line to the paragraph.
c
      if (line.ne.' ') call add_line (trim(line),current_par)
c
c ****** Deallocate storage.
c
      deallocate (done)
c
      return
      end
c#######################################################################
      subroutine defarg (group,keyword,default,name)
c
c-----------------------------------------------------------------------
c
c ****** Define the syntax for a command line argument item.
c
c-----------------------------------------------------------------------
c
c ****** GROUP is the syntax group index;
c ****** KEYWORD is the keyword;
c ****** DEFAULT is the default value of the argument;
c ****** NAME is the name of the argument (for use in
c ****** constructing the usage line).
c
c-----------------------------------------------------------------------
c
      use syntax
      use parse_args
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: group
      character(*) :: keyword
      character(*) :: default
      character(*) :: name
c
c-----------------------------------------------------------------------
c
      integer :: i,n,ix
      character(512), dimension(:), allocatable :: names
c
c-----------------------------------------------------------------------
c
      integer, external :: number_of_names
      integer, external :: match_keyword
c
c-----------------------------------------------------------------------
c
c ****** Check that the group index is valid.
c
      if (group.lt.0.or.group.gt.ngroups) then
        write (*,*)
        write (*,*) '### ERROR in DEFARG:'
        write (*,*) '### An invalid group index was specified.'
        write (*,*) 'Group index = ',group
        write (*,*) 'Keyword: ',trim(keyword)
        if (name.ne.' ') write (*,*) 'Name: ',trim(name)
        write (*,*)
        write (*,*) '### This indicates a programming error'//
     &              ' in the syntax definition and use.'
        call exit (1)
      end if
c
c ****** Check for a null keyword.
c
      if (keyword.eq.' ') then
        write (*,*)
        write (*,*) '### ERROR in DEFARG:'
        write (*,*) '### The keyword is null.'
        write (*,*) 'Group index = ',group
        if (name.ne.' ') write (*,*) 'Name: ',trim(name)
        write (*,*)
        write (*,*) '### This indicates a programming error'//
     &              ' in the syntax definition and use.'
        call exit (1)
      end if
c
c ****** Check if this defines a group-set.
c
      if (group.eq.GROUP_K_ONE_ONLY.or.
     &    group.eq.GROUP_K_ONE_OR_NONE) then
c
c ****** Increment the group-set counter.
c
        if (ngroup_sets.ge.mxgroup_sets) then
          write (*,*)
          write (*,*) '### ERROR in DEFARG:'
          write (*,*) '### Exceeded the number of allowed group-sets.'
          write (*,*) 'Maximum number of group-sets = ',mxgroup_sets
          write (*,*) 'Group index = ',group
          write (*,*) 'Keyword: ',trim(keyword)
          write (*,*)
          write (*,*) '### This indicates a programming error'//
     &                ' in the syntax definition and use.'
          call exit (1)
        end if
c
        ngroup_sets=ngroup_sets+1
c
c ****** Set the group type.
c
        group_sets(ngroup_sets)%group_type=group
c
c ****** Set the appropriate group-set flags.
c
        if (group.eq.GROUP_K_ONE_ONLY) then
          group_sets(ngroup_sets)%required=.true.
          group_sets(ngroup_sets)%only_one=.true.
        else if (group.eq.GROUP_K_ONE_OR_NONE) then
          group_sets(ngroup_sets)%required=.false.
          group_sets(ngroup_sets)%only_one=.true.
        else
          group_sets(ngroup_sets)%required=.false.
          group_sets(ngroup_sets)%only_one=.false.
        end if
c
c ****** Mark the members in the group-set.
c
        n=number_of_names(keyword)
        allocate (names(n))
        allocate (group_sets(ngroup_sets)%members(n))
        group_sets(ngroup_sets)%n_members=n
        call load_list_of_names (keyword,n,names)
        do i=1,n
          ix=match_keyword(names(i))
          if (ix.eq.0) then
            write (*,*)
            write (*,*) '### ERROR in DEFARG:'
            write (*,*) '### Error while defining a group set.'
            write (*,*) 'Could not match a group-set keyword to'//
     &                  ' an existing keyword:'
            write (*,*) 'Group-set keyword: ',trim(keyword)
            write (*,*) 'Keyword not matched: ',trim(names(i))
            write (*,*)
            write (*,*) '### This indicates a programming error'//
     &                  ' in the syntax definition and use.'
            call exit (1)
          end if
          group_sets(ngroup_sets)%members(i)=ix
        enddo
        deallocate (names)
c
        return
c
      end if
c
c ****** Increment the argument counter.
c
      if (nargs.ge.mxarg) then
        write (*,*)
        write (*,*) '### ERROR in DEFARG:'
        write (*,*) '### Exceeded the number of allowed arguments.'
        write (*,*) 'Maximum number of arguments = ',mxarg
        write (*,*) 'Group index = ',group
        write (*,*) 'Keyword: ',trim(keyword)
        if (name.ne.' ') write (*,*) 'Name: ',trim(name)
        write (*,*)
        write (*,*) '### This indicates a programming error'//
     &              ' in the syntax definition and use.'
        call exit (1)
      end if
c
      nargs=nargs+1
c
c ****** Store the group index and keyword.
c
c ****** For group GROUP_A (single arguments), the name of the
c ****** argument is passed as the "keyword".
c
      args(nargs)%group=group
      call put_str (trim(keyword),args(nargs)%keyword)
c
c ****** Initialize the flag that indicates whether an argument
c ****** has been set.
c
      args(nargs)%set=.false.
c
c ****** If a default argument was supplied, the argument
c ****** does not have to be set.  Use DEFAULT=' ' to
c ****** indicate that an argument is required.
c
c ****** If a default argument has been supplied, store it in
c ****** ARGS(nargs)%VALUE.  If there is no default,
c ****** set ARGS(nargs)%VALUE to an empty string.
c
c ****** Since group GROUP_K doesn't have an argument,
c ****** DEFAULT is ignored for this group.
c
      if (group.eq.GROUP_K) then
        args(nargs)%required=.false.
        call put_str (' ',args(nargs)%value)
      else
        if (default.eq.' ') then
          args(nargs)%required=.true.
          call put_str (' ',args(nargs)%value)
        else
          args(nargs)%required=.false.
          call put_str (trim(default),args(nargs)%value)
        end if
      end if
c
c ****** Store the argument name.  For groups GROUP_K (keywords)
c ****** and GROUP_A (single arguments), there is no argument name,
c ****** so NAME is ignored.
c
      if (group.eq.GROUP_K.or.group.eq.GROUP_A) then
        call put_str (' ',args(nargs)%name)
      else
        call put_str (trim(name),args(nargs)%name)
      end if
c
      return
      end
c#######################################################################
      subroutine fetcharg (keyword,set,arg)
c
c-----------------------------------------------------------------------
c
c ****** Fetch the value of the argument corresponding to
c ****** keyword KEYWORD.
c
c-----------------------------------------------------------------------
c
c ****** If KEYWORD is a keyword-type argument (GROUP_K), return
c ****** its setting through variable SET.  The variable ARG should
c ****** be ignored for this type of keyword.
c
c ****** For keywords with arguments (GROUP_A, GROUP_KA, and
c ****** GROUP_KAA), return the value of the arguments in ARG,
c ****** and return SET=.true. if they were set via the command line;
c ****** otherwise, return SET=.false..
c
c-----------------------------------------------------------------------
c
      use parse_args
      use get_str_interface
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(*) :: keyword
      logical :: set
      character(*) :: arg
c
c-----------------------------------------------------------------------
c
      integer :: i
c
c-----------------------------------------------------------------------
c
      do i=nargs,1,-1
        if (keyword.eq.get_str(args(i)%keyword)) go to 100
      enddo
c
      write (*,*)
      write (*,*) '### ERROR in FETCHARG:'
      write (*,*) '### The requested keyword could not be matched.'
      write (*,*) 'Keyword = ',trim(keyword)
      write (*,*)
      write (*,*) '### This indicates a programming error'//
     &            ' in the syntax definition and use.'
      call exit (1)
c
  100 continue
c
      set=args(i)%set
      arg=get_str(args(i)%value)
c
      return
      end
c#######################################################################
      function matchkw (group,keyword)
c
c-----------------------------------------------------------------------
c
c ****** Match keyword KEYWORD against the list of keywords in
c ****** group GROUP.
c
c ****** If found, set the function value to the corresponding
c ****** argument number.  Otherwise, return MATCHKW=0.
c
c-----------------------------------------------------------------------
c
      use parse_args
      use get_str_interface
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: group
      character(*) :: keyword
      integer :: matchkw
c
c-----------------------------------------------------------------------
c
      integer :: i
c
c-----------------------------------------------------------------------
c
      matchkw=0
c
      do i=nargs,1,-1
        if (group.eq.args(i)%group) then
          if (keyword.eq.get_str(args(i)%keyword)) then
            matchkw=i
            return
          end if
        end if
      enddo
c
      return
      end
c#######################################################################
      function nextarg (group)
c
c-----------------------------------------------------------------------
c
c ****** Find the position of the next argument in group GROUP
c ****** that has not been set.
c
c-----------------------------------------------------------------------
c
c ****** If an empty slot is found, set the function value
c ****** to the corresponding argument number.
c
c ****** Otherwise, return NXTARG=0.
c
c-----------------------------------------------------------------------
c
      use parse_args
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: group
      integer :: nextarg
c
c-----------------------------------------------------------------------
c
      integer :: i
c
c-----------------------------------------------------------------------
c
      nextarg=0
c
      do i=1,nargs
        if (group.eq.args(i)%group) then
          if (.not.args(i)%set) then
            nextarg=i
            return
          end if
        end if
      enddo
c
      return
      end
c#######################################################################
      subroutine nargs_specified (n)
c
c-----------------------------------------------------------------------
c
c ****** Return the number of arguments specified on the command
c ****** line.
c
c-----------------------------------------------------------------------
c
      use parse_args
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: n
c
c-----------------------------------------------------------------------
c
      n=nargs_spec
c
      return
      end
c#######################################################################
      function match_keyword (keyword)
c
c-----------------------------------------------------------------------
c
c ****** Match KEYWORD to the current list of defined keywords.
c
c ****** A successful match returns the index of the first matched
c ****** entry in TABLE; an unsuccessful match returns 0.
c
c-----------------------------------------------------------------------
c
      use parse_args
      use get_str_interface
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(*), intent(in) :: keyword
      integer :: match_keyword
c
c-----------------------------------------------------------------------
c
      integer :: i
c
c-----------------------------------------------------------------------
c
      match_keyword=0
c
      do i=1,nargs
        if (keyword.eq.get_str(args(i)%keyword)) then
          match_keyword=i
          return
        end if
      end do
c
      return
      end
c#######################################################################
      subroutine new_par (par)
c
c-----------------------------------------------------------------------
c
c ****** Initialize paragraph PAR.
c
c-----------------------------------------------------------------------
c
      use paragraph_def
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(paragraph), pointer :: par
c
c-----------------------------------------------------------------------
c
      allocate (par)
      nullify (par%line%c)
      nullify (par%next)
c
      return
      end
c#######################################################################
      subroutine delete_par (par)
c
c-----------------------------------------------------------------------
c
c ****** Delete paragraph PAR and deallocate its storage and that
c ****** of its linked lists.
c
c-----------------------------------------------------------------------
c
      use paragraph_def
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(paragraph), pointer :: par
c
c-----------------------------------------------------------------------
c
      type(paragraph), pointer :: current_par,previous_par
c
c-----------------------------------------------------------------------
c
      current_par=>par
c
      do
c
c ****** Deallocate the line buffer.
c
        call delete_str (current_par%line)
c
c ****** Set the pointer to the next line (if it has been defined).
c
        if (.not.associated(current_par%next)) exit
        previous_par=>current_par
        current_par=>current_par%next
        deallocate (previous_par)
c
      enddo
c
      deallocate (current_par)
c
      return
      end
c#######################################################################
      subroutine add_line (line,par)
c
c-----------------------------------------------------------------------
c
c ****** Add LINE to paragraph PAR.
c
c ****** On exit from this routine, PAR points to a new line,
c ****** and can be used to store the next line of text.
c
c-----------------------------------------------------------------------
c
      use paragraph_def
      use new_par_interface
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(*) :: line
      type(paragraph), pointer :: par
c
c-----------------------------------------------------------------------
c
c ****** Store LINE into the string buffer for the current line.
c
      call put_str (line,par%line)
c
c ****** Allocate a pointer to the next line.
c
      call new_par (par%next)
c
c ****** Set PAR to point to the next line.
c
      par=>par%next
c
      return
      end
c#######################################################################
      subroutine print_par (par)
c
c-----------------------------------------------------------------------
c
c ****** Print all lines of paragraph PAR to STDOUT.
c
c-----------------------------------------------------------------------
c
      use paragraph_def
      use get_str_interface
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(paragraph), pointer :: par
c
c-----------------------------------------------------------------------
c
      type(paragraph), pointer :: current_par
c
c-----------------------------------------------------------------------
c
      current_par=>par
c
      do
c
c ****** Print the line if it has been defined.
c
        if (associated(current_par%line%c)) then
          write (*,*) trim(get_str(current_par%line))
        end if
c
c ****** Set the pointer to the next line (if it has been defined).
c
        if (.not.associated(current_par%next)) exit
        current_par=>current_par%next
c
      enddo
c
      return
      end
c#######################################################################
      subroutine put_str (cval,str)
c
c-----------------------------------------------------------------------
c
c ****** Store character variable CVAL into string STR.
c ****** This routine allocates storage for the string.
c
c-----------------------------------------------------------------------
c
      use string_def
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(*) :: cval
      type(string) :: str
c
c-----------------------------------------------------------------------
c
      integer :: l,i
c
c-----------------------------------------------------------------------
c
      l=len(cval)
c
      allocate (str%c(l))
c
      do i=1,l
        str%c(i)=cval(i:i)
      enddo
c
      return
      end
c#######################################################################
      function get_str (str)
c
c-----------------------------------------------------------------------
c
c ****** Return the value of string STR as the function value
c ****** (as an assumed-length character variable).
c
c-----------------------------------------------------------------------
c
      use string_def
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(string) :: str
      character(size(str%c)) :: get_str
c
c-----------------------------------------------------------------------
c
      integer :: i
c
c-----------------------------------------------------------------------
c
      do i=1,size(str%c)
        get_str(i:i)=str%c(i)
      enddo
c
      return
      end
c#######################################################################
      subroutine delete_str (str)
c
c-----------------------------------------------------------------------
c
c ****** Delete the storage for string STR.
c
c-----------------------------------------------------------------------
c
      use string_def
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(string) :: str
c
c-----------------------------------------------------------------------
c
      if (associated(str%c)) then
        deallocate (str%c)
      end if
      nullify (str%c)
c
      return
      end
c#######################################################################
      function intval (avalue,name)
c
c-----------------------------------------------------------------------
c
c ****** Get the value of the integer in character variable AVALUE.
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(*) :: avalue
      character(*) :: name
      integer :: intval
c
c-----------------------------------------------------------------------
c
      logical, external :: ifint
c
c-----------------------------------------------------------------------
c
      integer :: ivalue
c
c-----------------------------------------------------------------------
c
      if (.not.ifint(trim(avalue),ivalue)) then
        write (*,*)
        write (*,*) '### ERROR in INTVAL:'
        write (*,*) '### Could not interpret an integer '//
     &              'while setting: ',trim(name)
        write (*,*) 'Invalid format: ',trim(avalue)
        call exit (1)
      end if
c
      intval=ivalue
c
      return
      end
c#######################################################################
      function fpval (avalue,name)
c
c-----------------------------------------------------------------------
c
c ****** Get the value of the floating point number in character
c ****** variable AVALUE.
c
c-----------------------------------------------------------------------
c
      use number_types
c
c-----------------------------------------------------------------------
c
      character(*) :: avalue
      character(*) :: name
      real(r_typ) :: fpval
c
c-----------------------------------------------------------------------
c
      logical, external :: iffp
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: value
c
c-----------------------------------------------------------------------
c
      if (.not.iffp(trim(avalue),value)) then
        write (*,*)
        write (*,*) '### ERROR in FPVAL:'
        write (*,*) '### Could not interpret a floating point '//
     &              'number while setting: ',trim(name)
        write (*,*) 'Invalid format: ',trim(avalue)
        call exit (1)
      end if
c
      fpval=value
c
      return
      end
c#######################################################################
      function iffp (alpha,value)
c
c-----------------------------------------------------------------------
c
c ****** Determine if ALPHA represents a floating point number;
c ****** if so, return its value in VALUE.
c
c-----------------------------------------------------------------------
c
c ****** Set IFFP=.TRUE. if ALPHA contains an alphanumeric
c ****** string with the following format:
c
c       ALPHA = '[A][B...B][.][B...B][e[A]B[B...B]]',
c
c ****** where A represents a + or - sign, and B represents a digit
c ****** between 0 and 9, inclusive.
c ****** The exponent may be denoted by a lower or upper case e.
c ****** The mantissa must have at least one digit, and the
c ****** the exponent, if present, must have between 1 and 3 digits.
c ****** Otherwise, set IFFP=.FALSE.
c
c-----------------------------------------------------------------------
c
      use number_types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(*) :: alpha
      real(r_typ) :: value
      logical :: iffp
c
c-----------------------------------------------------------------------
c
      integer :: nmant,nexp,k,i,ke
      logical :: ifpoint,ifexp
      character(7) :: fmt
c
c-----------------------------------------------------------------------
c
      iffp=.false.
      ifpoint=.false.
      ifexp=.false.
      nmant=0
      nexp=0
c
      do k=1,len_trim(alpha)
        i=iachar(alpha(k:k))
c
c ****** Check for a sign in the first position.
c
        if (k.eq.1.and.(i.eq.43.or.i.eq.45)) cycle
c
c ****** Check for a digit.
c
        if (i.ge.48.and.i.le.57) then
c
c ****** Count digits in mantissa and exponent.
c
        if (ifexp) then
          nexp=nexp+1
          else
            nmant=nmant+1
          end if
          cycle
c
        end if
c
c ****** Check for a decimal point.
c
        if (.not.ifpoint.and.i.eq.46) then
c
c ****** Check that we are in the mantissa.
c
          if (.not.ifexp) then
            ifpoint=.true.
            cycle
          end if
c
        end if
c
c ****** Check for an exponent.
c
        if (.not.ifexp.and.(i.eq.101.or.i.eq.69)) then
          ifexp=.true.
          ke=k
          cycle
        end if
c
c ****** Check for an exponent sign.
c
        if (ifexp.and.k.eq.ke+1.and.(i.eq.43.or.i.eq.45)) cycle
c
c ****** Failed check: fall through here.
c
        iffp=.false.
c
        return
c
      enddo
c
c ****** Final check of validity: check number of digits in
c ****** the mantissa and exponent.
c
      if (nmant.ge.1) iffp=.true.
      if (ifexp.and.(nexp.lt.1.or.nexp.gt.3)) iffp=.false.
c
c ****** Obtain its numeric value.
c
      fmt='(f  .0)'
      write (fmt(3:4),'(i2.2)') len_trim(alpha)
c
      if (iffp) read (alpha,fmt) value
c
      return
      end
c#######################################################################
      function ifint (alpha,ivalue)
c
c-----------------------------------------------------------------------
c
c ****** If ALPHA represents an integer, return IFINT=.true., and
c ****** put its value into IVALUE.
c
c ****** Otherwise, return IFINT=.false..
c
c-----------------------------------------------------------------------
c
c ****** A valid integer has the format:
c
c          ALPHA = '[A]B[B...B]',
c
c ****** where A represents a + or - sign, and B represents a digit
c ****** between 0 and 9, inclusive.
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(*) :: alpha
      integer :: ivalue
      logical :: ifint
c
c-----------------------------------------------------------------------
c
      integer :: k,i
      character(5) :: fmt
c
c-----------------------------------------------------------------------
c
      ifint=.false.
c
      do k=1,len_trim(alpha)
c
        i=iachar(alpha(k:k))
c
c ****** Check for a sign in the first position.
c
        if (k.eq.1.and.(i.eq.43.or.i.eq.45)) cycle
c
c ****** Check for a digit.
c
        if (i.ge.48.and.i.le.57) then
          ifint=.true.
          cycle
        end if
c
c ****** Failed check: fall through here.
c
        ifint=.false.
c
        return
c
      enddo
c
c ****** Obtain its numeric value.
c
      fmt='(i  )'
      write (fmt(3:4),'(i2.2)') len_trim(alpha)
c
      if (ifint) read (alpha,fmt) ivalue
c
      return
      end
c#######################################################################
      subroutine load_list_of_reals (s,label,n,f)
c
c-----------------------------------------------------------------------
c
c ****** Read N real values from character string S into
c ****** array F(N). The values in S may be either space or
c ****** comma separated.
c
c-----------------------------------------------------------------------
c
      use number_types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(*) :: s
      character(*) :: label
      integer :: n
      real(r_typ), dimension(n) :: f
      intent(in) :: s,label,n
      intent(out) :: f
c
c-----------------------------------------------------------------------
c
      integer :: i,i0,i1
      character :: delimiter
      character(512) :: list
c
c-----------------------------------------------------------------------
c
      real(r_typ), external :: fpval
c
c-----------------------------------------------------------------------
c
c ****** Make a local copy of the string (removing leading spaces).
c
      list=adjustl(s)
c
c ****** If any commas are present, use a comma as the delimiter.
c ****** Otherwise, one or more spaces is used as a delimiter.
c ****** In this case, compress multiple spaces into a single space.
c
      if (index(list,',').ne.0) then
        delimiter=','
      else
        delimiter=' '
        call delete_repeated_char (list,' ')
      end if
c
c ****** Read the list of N numbers sequentially into F.
c
      i0=1
      do i=1,n-1
        i1=scan(list(i0:),delimiter)+i0-2
        f(i)=fpval(adjustl(list(i0:i1)),label)
        i0=i1+2
      enddo
      f(n)=fpval(adjustl(list(i0:)),label)
c
      return
      end
c#######################################################################
      subroutine load_list_of_ints (s,label,n,j)
c
c-----------------------------------------------------------------------
c
c ****** Read N integer values from character string S into
c ****** array J(N).  The values in S may be either space or
c ****** comma separated.
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(*) :: s
      character(*) :: label
      integer :: n
      integer, dimension(n) :: j
      intent(in) :: s,label,n
      intent(out) :: j
c
c-----------------------------------------------------------------------
c
      integer :: i,i0,i1
      character :: delimiter
      character(512) :: list
c
c-----------------------------------------------------------------------
c
      integer, external :: intval
c
c-----------------------------------------------------------------------
c
c ****** Make a local copy of the string (removing leading spaces).
c
      list=adjustl(s)
c
c ****** If any commas are present, use a comma as the delimiter.
c ****** Otherwise, one or more spaces is used as a delimiter.
c ****** In this case, compress multiple spaces into a single space.
c
      if (index(list,',').ne.0) then
        delimiter=','
      else
        delimiter=' '
        call delete_repeated_char (list,' ')
      end if
c
c ****** Read the list of N numbers sequentially into J.
c
      i0=1
      do i=1,n-1
        i1=scan(list(i0:),delimiter)+i0-2
        j(i)=intval(adjustl(list(i0:i1)),label)
        i0=i1+2
      enddo
      j(n)=intval(adjustl(list(i0:)),label)
c
      return
      end
c#######################################################################
      function number_of_names (s)
c
c-----------------------------------------------------------------------
c
c ****** Return the number of names in string S.  The values
c ****** can be separated by either spaces or commas.
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(*) :: s
      integer :: number_of_names
      intent(in) :: s
c
c-----------------------------------------------------------------------
c
      integer :: l,n,i0,i1
      character :: delimiter
      character(len(s)) :: list
c
c-----------------------------------------------------------------------
c
c ****** Make a local copy of the string (removing leading spaces).
c
      list=adjustl(s)
c
c ****** If any commas are present, use a comma as the delimiter.
c ****** Otherwise, one or more spaces is used as a delimiter.
c ****** In this case, compress multiple spaces into a single space.
c
      if (index(list,',').ne.0) then
        delimiter=','
      else
        delimiter=' '
        call delete_repeated_char (list,' ')
      end if
c
c ****** Find the number of names N.
c
      l=len_trim(list)
c
      n=0
c
      i0=1
      do
        i1=scan(list(i0:l),delimiter)
        if (i1.eq.0) then
          if (.not.(delimiter.eq.' '.and.l.lt.i0)) then
            n=n+1
          end if
          exit
        end if
        i1=i1+i0-2
        i0=i1+2
        n=n+1
      enddo
c
      number_of_names=n
c
      return
      end
c#######################################################################
      subroutine load_list_of_names (s,n,names)
c
c-----------------------------------------------------------------------
c
c ****** Read up to N names from character string S into
c ****** array NAMES(N).  The values can be separated by either
c ****** spaces or commas.
c
c-----------------------------------------------------------------------
c
c ****** This routine is designed to be used in conjunction with
c ****** routine NUMBER_OF_NAMES.  A typical use would be:
c
c          character(32) :: s
c          character(32), dimension(:), allocatable :: names
c          integer :: n
c
c          n=number_of_names(s)
c          allocate (names(n))
c          call load_list_of_names (s,n,names)
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(*) :: s
      integer :: n
      character(*), dimension(n) :: names
      intent(in) :: s,n
      intent(out) :: names
c
c-----------------------------------------------------------------------
c
      integer :: l,i,i0,i1
      character :: delimiter
      character(len(s)) :: list
c
c-----------------------------------------------------------------------
c
c ****** Make a local copy of the string (removing leading spaces).
c
      list=adjustl(s)
c
c ****** If any commas are present, use a comma as the delimiter.
c ****** Otherwise, one or more spaces is used as a delimiter.
c ****** In this case, compress multiple spaces into a single space.
c
      if (index(list,',').ne.0) then
        delimiter=','
      else
        delimiter=' '
        call delete_repeated_char (list,' ')
      end if
c
c ****** Read the names.
c
      l=len_trim(list)
c
      i=0
c
      i0=1
      do
        i1=scan(list(i0:l),delimiter)
        if (i1.eq.0) then
          if (.not.(delimiter.eq.' '.and.l.lt.i0)) then
            i=i+1
            if (i.gt.n) exit
            names(i)=adjustl(list(i0:l))
          end if
          exit
        end if
        i=i+1
        if (i.gt.n) exit
        i1=i1+i0-2
        names(i)=adjustl(list(i0:i1))
        i0=i1+2
      enddo
c
      return
      end
c#######################################################################
      subroutine delete_repeated_char (s,c)
c
c-----------------------------------------------------------------------
c
c ****** Transform repeated adjoining occurrences of character C
c ****** in string S into single occurrences of C.
c
c ****** The string S is overwritten by the modified string.
c
c ****** Trailing blanks in S are ignored.
c
c-----------------------------------------------------------------------
c
c ****** For example, suppose this routine is called with C='d' and
c ****** S='abcdddeefdhdd'.  On return, S will have the value
c ****** 'abcdeefdhd'.
c
c-----------------------------------------------------------------------
c
c ****** This routine uses the FORTRAN90 intrinsic SCAN.
c
c-----------------------------------------------------------------------
c
      character(*) :: s
      character :: c
      intent(in) :: c
      intent(inout) :: s
c
c-----------------------------------------------------------------------
c
      integer :: i,i0
c
c-----------------------------------------------------------------------
c
      i0=1
      do
        i=scan(trim(s(i0:)),c)
        if (i.eq.0) exit
        i0=i0+i
        do
          if (s(i0:i0).ne.c) exit
          s(i0:)=s(i0+1:)
        enddo
      enddo
c
      return
      end
