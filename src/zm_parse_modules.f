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
c#######################################################################
      module syntax
c
c-----------------------------------------------------------------------
c ****** Group definitions for parsing command-line arguments.
c-----------------------------------------------------------------------
c
c        GROUP 1: <kw>
c        GROUP 2: <arg>
c        GROUP 3: <kw> <arg>
c        GROUP 4: <kw> <arg> <arg>
c        GROUP 5: GROUP SET (only one must be specified)
c        GROUP 6: GROUP SET (only one or none must be specified)
c
      integer, parameter :: ngroups=6
c
      integer, parameter :: GROUP_K            =1
      integer, parameter :: GROUP_A            =2
      integer, parameter :: GROUP_KA           =3
      integer, parameter :: GROUP_KAA          =4
      integer, parameter :: GROUP_K_ONE_ONLY   =5
      integer, parameter :: GROUP_K_ONE_OR_NONE=6
c
      end module
c#######################################################################
      module string_def
c
c-----------------------------------------------------------------------
c ****** Define a structure to hold a string.
c-----------------------------------------------------------------------
c
      implicit none
c
      type :: string
        character, dimension(:), pointer :: c
      end type
c
      end module
c#######################################################################
      module paragraph_def
c
      use string_def
c
      implicit none
c
c-----------------------------------------------------------------------
c ****** Define a structure for a linked list of lines
c ****** (i.e., a paragraph).
c-----------------------------------------------------------------------
c
      type :: paragraph
        type(string) :: line
        type(paragraph), pointer :: next
      end type
c
c-----------------------------------------------------------------------
c ****** Define a structure to hold a list of paragraphs.
c-----------------------------------------------------------------------
c
      type :: parlist
        type(paragraph), pointer :: par
      end type
c
      end module
c#######################################################################
      module lcase_interface
      interface
        function lcase (s)
        character(*) :: s
        character(len(s)) :: lcase
        end function
      end interface
      end module
c#######################################################################
      module ucase_interface
      interface
        function ucase (s)
        character(*) :: s
        character(len(s)) :: ucase
        end function
      end interface
      end module
c#######################################################################
      module new_par_interface
      interface
        subroutine new_par (par)
        use paragraph_def
        implicit none
        type(paragraph), pointer :: par
        end subroutine
      end interface
      end module
c#######################################################################
      module delete_par_interface
      interface
        subroutine delete_par (par)
        use paragraph_def
        implicit none
        type(paragraph), pointer :: par
        end subroutine
      end interface
      end module
c#######################################################################
      module add_line_interface
      interface
        subroutine add_line (line,par)
        use paragraph_def
        implicit none
        character(*) :: line
        type(paragraph), pointer :: par
        end subroutine
      end interface
      end module
c#######################################################################
      module print_par_interface
      interface
        subroutine print_par (par)
        use paragraph_def
        implicit none
        type(paragraph), pointer :: par
        end subroutine
      end interface
      end module
c#######################################################################
      module get_str_interface
      interface
        function get_str (str)
        use string_def
        implicit none
        type(string) :: str
        character(size(str%c)) :: get_str
        end function
      end interface
      end module
c#######################################################################
      module get_usage_line_interface
      interface
        subroutine get_usage_line (usage)
        use paragraph_def
        implicit none
        type(paragraph), pointer :: usage
        end subroutine
      end interface
      end module
