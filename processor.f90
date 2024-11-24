module processor

    !! This module contains input and output subroutines
    !!
    !! ** Note: It is not necessary to alter in this module **

    use fedata
    implicit none

    public :: input, output, stopwatch, plotmatlabdef, plotmatlabeval, plotmatlabevec, plotmatlabeig
    private

    ! Counters for output files
    integer :: imatlab = 0


contains

    subroutine input

        !!  This subroutine reads in the input in ansys format.

        integer :: k, e
        integer :: nncheck, necheck
        integer :: net
        integer :: ivalue, ivalue2
        integer :: currentet, currentmp, currentr
        character(len=12), dimension(10) :: ename
        character(len=12) :: command, cvalue, cvalue2
        real(wp) :: rvalue, rvalue2
        logical :: fnexist
        integer, dimension(:,:), allocatable :: eix

        ! Reset counters
        nn = 0
        nncheck = 0
        net = 0
        ne = 0
        necheck = 0
        nm = 0
        nb = 0
        np = 0

        ! Note: From Fortran 2003+ command line arguments are accessible, but not in Fortran90
        write (*, *)
        write (*, '("File(s) in current directory: ")')
        call system('dir') ! Windows
!        call system('ls' ) ! Linux
        write (*, *)
        write (*, '("Enter fem input file (<filename> or enter for previous <filename>)")')
        write (*, '("Filename?")')
        read (*, '(a50)') filename
        if (filename == ' ') then
            inquire (file = '.fem_filename', exist = fnexist)
            if (.not. fnexist) then
                write (*, *)
                write (*, '("ERROR: previous filename does not exist")')
                stop
            end if
            open (10, file = '.fem_filename')
            read (10, '(a50)') filename
            close (10)
        end if
        inquire (file = trim(filename), exist = fnexist)
        if (.not. fnexist) then
            write (*, *)
            write (*, '("ERROR: file ", a, " does not exist")') trim(filename)
            stop
        end if

        open (10, file = '.fem_filename')
        write (10, '(a)') trim(filename)
        close (10)

        write (*, *)
        write (*, '("Reading input file ", a)') trim(filename)

        open (10, file = trim(filename))

        ! Pass 1: Scan through file in order to find out how many
        ! elements, nodes, etc that are defined

        ! Read in /PREP7 data
        do
            read (10, '(a)', end=100) command
            if (command == '/PREP7') exit
            if (command == '/' .or. command == ' ' .or. command == '!' ) cycle
        end do
        do
            read (10, '(a1)', end = 200) command
            if (command == '/' .or. command == ' ') cycle

            backspace (10)
            read (10, *) command
            select case (command)
            case ('FINISH', 'finish')
                exit
            case ('N', 'n')
                nn = nn + 1
                backspace (10)
                read (10, *) command, ivalue
                if (ivalue > nncheck) nncheck = ivalue
            case('ET', 'et')
                net = net + 1
            case('EN', 'en')
                ne = ne + 1
                backspace (10)
                read (10, *) command, ivalue
                if (ivalue > necheck) necheck = ivalue
            case('MP', 'mp')
                backspace (10)
                read (10, *) command, cvalue, ivalue
                if (ivalue > nm) nm = ivalue
            case('D', 'd')
                nb = nb + 1
            case('F', 'f')
                np = np + 1
            case('SFE', 'sfe')
                np = np + 1
            end select
        end do

        if (nn == 0) then
            write (*, *) 'ERROR: No nodes defined'
            stop
        else if (nn /= nncheck) then
            write (*, *) 'ERROR: Node number(s) skipped'
            stop
        else if (net == 0) then
            write (*, *) 'ERROR: No element types defined'
            stop
        else if (ne == 0) then
            write (*, *) 'ERROR: No elements defined'
            stop
        else if (ne /= necheck) then
            write (*, *) 'ERROR: Element number(s) skipped'
            stop
        else if (nm == 0) then
            write (*, *) 'ERROR: No material properties defined'
            stop
        else if (nb == 0) then
            write (*, *) 'ERROR: No supports defined'
            stop
        else if (np == 0) then
            write (*, *) 'WARNING: No loads defined'
        end if
        close (10)

        allocate (x(nn, 2))
        allocate (element(ne))
        allocate (eix(ne,4))
        allocate (mprop(nm))
        allocate (bound(nb, 3), loads(np, 4))

        ! Make thickness equal to one as default
        mprop%thk = 1

        ! And reset all other mprop parameters
        mprop%young = 0
        mprop%nu = 0
        mprop%dens = 0
        mprop%youngy = 0
        mprop%shear = 0
        do e = 1, ne
            element(e)%ix = 0
        end do

        open (10, file = trim(filename))

        nb = 0
        np = 0
        currentet = 0
        currentmp = 0
        currentr = 0
        accel=0.0

        ! Pass 2: Read model definition
        do
            read (10, '(a1)', end=200) command
            if (command == '/' .or. command == ' ' .or. command == '!' ) cycle

            backspace (10)
            read (10, *) command
            select case (command)
            case ('FINISH', 'finish')
                exit
            case ('N', 'n')
                backspace (10)
                read (10, *) command, ivalue, x(ivalue, 1), x(ivalue, 2), rvalue

            case ('ET', 'et')
                backspace (10)
                read (10, *) command, ivalue, cvalue
                select case (cvalue)
                case ('LINK1',  'link1')
                    ename(ivalue) = 'link1'
                case ('PLANE42', 'plane42')
                    ename(ivalue) = 'plane42'
                case ('PLANE42RECT', 'plane42rect')
                    ename(ivalue) = 'plane42'
                case default
                    write (*, *) 'ERROR: Undefined element: ', trim(cvalue)
                    stop
                end select
                currentet = ivalue
            case ('TYPE', 'type')
                backspace (10)
                read (10, *) command, ivalue
                currentet = ivalue
            case ('EN', 'en')
                if (currentet == 0) then
                    write (*, *) 'ERROR: No previous element type pointer defined'
                    stop
                end if
                if (currentmp == 0) then
                    write (*, *) 'ERROR: No previous material property pointer defined'
                    stop
                end if
                backspace (10)
                select case (ename(currentet))
                case ('link1')
                    read (10, *) command, ivalue, (element(ivalue)%ix(k), k = 1, 2)
                    element(ivalue)%mat = currentmp
                    element(ivalue)%id  = 1
                    element(ivalue)%numnode = 2
                case ('plane42')
                    read (10, *) command, ivalue, (element(ivalue)%ix(k), k = 1, 4)
                    element(ivalue)%mat = currentmp
                    element(ivalue)%id  = 2
                    element(ivalue)%numnode = 4
                case default
                    write(*, *) 'ERROR: Unknown element type'
                    stop
                end select
            case ('MP', 'mp')
                backspace (10)
                read (10, *) command, cvalue, ivalue, rvalue
                select case (cvalue)
                case ('EX', 'ex')
                    mprop(ivalue)%young = rvalue
                case ('PRXY', 'prxy')
                    mprop(ivalue)%nu = rvalue
                case ('DENS', 'dens')
                    mprop(ivalue)%dens = rvalue
                case ('EY', 'ey')
                    mprop(ivalue)%youngy = rvalue
                case ('GXY', 'gxy')
                    mprop(ivalue)%shear = rvalue
                case default
                    write (*, *) 'ERROR: Undefined material property: ', trim(cvalue)
                    stop
                end select
                currentmp = ivalue
            case ('MAT', 'mat')
                backspace (10)
                read (10, *) command, ivalue
                currentmp = ivalue
            case ('R', 'r')
                backspace (10)
                read (10, *) command, ivalue, rvalue
                select case (ename(currentet))
                case ('link1')
                    ! define truss area
                    mprop(ivalue)%area = rvalue
                case ('plane42')
                    ! define 4-noded quad element thickness
                    mprop(ivalue)%thk = rvalue
                case default
                    write (*, *) 'ERROR: Undefined real constant (r) card'
                    stop
                end select
                currentr = ivalue
            case ('REAL', 'real')
                ! TODO
                ! This command is only partially implemented -- it should not be used
                ! (it has no effect at the moment). The command is supposed to set
                ! the default value of "R", which is to be used for area or thickness
                ! in materials, that do not have this value set explicitly (via a
                ! call of "R" command before calling "MP" on some card).
                backspace (10)
                read (10, *) command, ivalue
                currentr = ivalue
            case ('D', 'd')
                backspace (10)
                nb = nb + 1
                read (10, *) command, ivalue, cvalue, rvalue
                bound(nb, 1) = ivalue
                select case (cvalue)
                case ('UX', 'ux')
                    bound(nb, 2) = 1
                case ('UY', 'uy')
                    bound(nb, 2) = 2
                case default
                    write (*, *) 'ERROR: Nodal dof not "ux" or "uy"'
                    stop
                end select
                bound(nb, 3) = rvalue
            case ('F', 'f')
                backspace (10)
                np = np + 1
                read (10, *) command, ivalue, cvalue, rvalue
                loads(np, 1) = 1
                loads(np, 2) = ivalue
                select case (cvalue)
                case ('FX', 'fx')
                    loads(np, 3) = 1
                case ('FY', 'fy')
                    loads(np, 3) = 2
                case default
                    write (*, *) 'ERROR: Nodal dof not "fx" or "fy"'
                    stop
                end select
                loads(np, 4) = rvalue
            case ('SFE', 'sfe')
                backspace (10)
                np = np + 1
                read (10, *) command, ivalue, ivalue2, cvalue, cvalue2, rvalue
                loads(np, 1) = 2
                loads(np, 2) = ivalue
                loads(np, 3) = ivalue2
                loads(np, 4) = rvalue
            case ('ACEL', 'acel')
                backspace (10)
                read (10, *) command, rvalue, rvalue2
                accel(1) = rvalue
                accel(2) = rvalue2
            case default
                write (*, *) 'ERROR: Unknown keyword: ', trim(command)
                stop
            end select
        end do

        close (10)

        ! Set default analysis type
        antype = 'STATIC'

        open (10, file = trim(filename))

        ! Read in /SOLU data (if it exists)
        do
          read (10, '(a)', end=400) command
          if (command == '/SOLU' .or. command == '/solu' ) exit
          if (command == '/' .or. command == ' ') cycle
        end do
        do
          read (10, '(a1)', end=300) command
          if (command == '/' .or. command == ' ') cycle

          backspace (10)
          read (10, *) command
          select case (command)
          case ('FINISH', 'finish')
                exit
          case ('ANTYPE', 'antype')
                backspace (10)
                read (10, *) command, cvalue
                select case (cvalue)
                case ('STATIC', 'static')
                    antype = 'static'
                case ('STATIC_NL', 'static_nl')
                    antype = 'static_nl'
                case ('MODAL', 'modal')
                    antype = 'modal'
                case ('ANGLE', 'angle')
                    antype = 'angle'
                case ('TRANS', 'trans')
                    antype = 'trans'
                case default
                    write (*, *) 'ERROR: Only static, static_nl, modal, angle and trans analyses are implemented'
                    stop
                end select
            end select
        end do

        400 continue

        return

        100 write (*, *)
        write (*, '("ERROR: no /PREP7 in input file")')
        stop

        200 write (*, *)
        write (*, '("ERROR: no /PREP7 FINISH in input file")')
        stop

        300 write (*, *)
        write (*, '("ERROR: no /SOLU FINISH in input file")')
        stop
    end subroutine input
!
!--------------------------------------------------------------------------------------------------
!
    subroutine output

        !!  This subroutine writes results into a text file

        integer :: i, e

        write (*, *)
        write (*, '("Writing fem output file to ", a)') trim(filename)//'.txt'
        open (10, file=trim(filename)//'.txt')

        write (10, '(" Node            Displacements                Applied/reaction forces")')
        write (10, '("number     1-direction     2-direction      1-direction     2-direction")')
        do i = 1, nn
            write (10, '(i6, 1x, f15.9, 1x, f15.9, 2x, f15.9, 1x, f15.9)') i, d(2*i-1), d(2*i), p(2*i-1), p(2*i)
        end do
        write (10, '("                                            ___________     ___________")')
        write (10, '("                              Summation ", f15.9, 1x, f15.9)') &
          sum(p(1:neqn-1:2)), sum(p(2:neqn:2))

        write (10, *)
        write (10, '("Element                    Element strain                                   Element Stress")')
        write (10, '("number      1-direction     2-direction     3-direction      1-direction     2-direction     3-direction")')
        do e = 1, ne
            select case( element(e)%id )
            case( 1 )
                write (10, '(1x, i6, 1x, f15.9, 34x, f15.9)') e, strain(e, 1), stress(e, 1)
            case( 2 )
                write (10, '(1x, i6, 3(1x, f15.9), 1x, 3(1x, f15.9))') e, (strain(e, i), i = 1, 3) ,(stress(e, i), i = 1, 3)
            end select
        end do

        close (10)
    end subroutine output
!
!--------------------------------------------------------------------------------------------------
!
    subroutine stopwatch(oper)

        !! This subroutine computes elapsed wallclock time
        !!
        !! Timing information is written in terminal window

        character(len=4), intent(in) :: oper
            !! Select which "button" to press on your Stopwatch:
            !!
            !! * 'STAR' or 'star' = reset and start the stopwatch
            !! * 'STOP' or 'stop' = print time spent since last 'star' operation (the stop watch is not reset)
        integer time_array_0(8), time_array_1(8)
        real(wp), save :: start_time, finish_time

        select case (oper)
        case ('STAR', 'star')
            call date_and_time(values=time_array_0)
            start_time = time_array_0 (5) * 3600 + time_array_0 (6) * 60 + time_array_0 (7) + 0.001 * time_array_0 (8)
        case ('STOP', 'stop')
            call date_and_time(values=time_array_1)
            finish_time = time_array_1 (5) * 3600 + time_array_1 (6) * 60 + time_array_1 (7) + 0.001 * time_array_1 (8)
            write (6, '(8x, 1a, 1f16.6)') 'elapsed wall clock time:', finish_time - start_time
        case default
            write (*, '("ERROR: in Processor/stopwatch")')
            stop
        end select
    end subroutine stopwatch

!
!--------------------------------------------------------------------------------------------------
!
    subroutine plotmatlabdef(title)

        ! Subroutine to plot the un/deformed structure using Matlab

        character(len=*), intent(in) :: title
        character(len=100) :: fdata, fscript
        integer i, j, e

        ! write datafile
        imatlab = imatlab + 1
        write (fdata, '(a,"_plotdeformed_data",i0.4,".m")') trim(filename), imatlab
        open (13, file=trim(fdata))

        ! write nodal coordinates
        write(13,'("X = [")')
        do i = 1,size(x,1)
            write (13,'(3(f15.9,1x))') (x(i,j), j=1,size(x,2))
        end do
        write(13,'("];")')
        write(13,'( )')

        ! write topology matrix
        write(13,'("IX = [")')
        do e = 1, ne
            write (13,'(20(i6,1x))') (element(e)%ix(i) ,i=1,size(element(e)%ix,1))
        end do
        write(13,'("];")')

        ! write deformation vector
        write(13,'("D = [")')
        do i = 1, neqn
            write (13,'(f15.9,1x)') (d(i) )
        end do
        write(13,'("];")')
        close(13)

        ! Create matlab script
        write(fscript, '(a,"_plotdeformed_",i0.4,".m")' ) trim(filename), imatlab
        open (13, file = trim(fscript))
        write(13,*) '% Plotting Un-Deformed and Deformed Structure'
        write(13,*) 'close all'
        write(13,*) 'clear all'
        write(13,*) fdata(1:len_trim(fdata)-2) // ';'
        write(13,*) '% Make plot'
        write(13,*) 'figure'
        write(13,*) 'set(gcf,',"'",'Name',"','", trim(title) ,"'",')'

        ! Element dependent code
        ! NOTE: not possible to mix element types !!!
        if (element(1)%id == 1) then
            write(13,*) 'subplot(2,1,2)'
            write(13,*) 'hold on'
            write(13,*) 'for e = 1:size(IX,1)'
            write(13,*) '   edof = [2*IX(e,1)-1 2*IX(e,1) 2*IX(e,2)-1 2*IX(e,2)];'
            write(13,*) '   xx = X(IX(e,1:2),1) + D(edof(1:2:4));'
            write(13,*) '   yy = X(IX(e,1:2),2) + D(edof(2:2:4));'
            write(13,*) '   plot(xx,yy,',"'",'b',"',","'",'LineWidth',"'",',1.5)'
            write(13,*) 'end'
            write(13,*) 'title(',"'",'Deformed',"'",')'
            write(13,*) 'axis equal'
            write(13,*) 'xaxes = get(gca,',"'",'xlim',"'",');'
            write(13,*) 'yaxes = get(gca,',"'",'ylim',"'",');'
            write(13,*) 'axis off;'
            write(13,*) 'subplot(2,1,1)'
            write(13,*) 'hold on'
            write(13,*) 'for e = 1:size(IX,1)'
            write(13,*) '   xx = X(IX(e,1:2),1);'
            write(13,*) '   yy = X(IX(e,1:2),2);'
            write(13,*) '   plot(xx,yy,',"'",'b',"',","'",'LineWidth',"',",'1.5)'
            write(13,*) 'end'
            write(13,*) 'title(',"'",'Undeformed',"'",')'
            write(13,*) 'axis([min(xaxes(1),min(X(:,1))) max(xaxes(2),max(X(:,1)))...'
            write(13,*) ' min(yaxes(1),min(X(:,2))) max(yaxes(2),max(X(:,2))) ]);'
            write(13,*) 'axis off;'
        else if (element(1)%id == 2) then
            write(13,*) 'subplot(2,1,2)'
            write(13,*) 'hold on'
            write(13,*) 'for e = 1:size(IX,1)'
            write(13,*) '    edof = [2*IX(e,1)-1 2*IX(e,1) 2*IX(e,2)-1 2*IX(e,2)...'
            write(13,*) '        2*IX(e,3)-1 2*IX(e,3) 2*IX(e,4)-1 2*IX(e,4)];'
            write(13,*) '    xx = X(IX(e,1:4),1) + D(edof(1:2:8));'
            write(13,*) '    yy = X(IX(e,1:4),2) + D(edof(2:2:8));'
            write(13,*) '    patch(xx,yy,[1 1 0]);'
            write(13,*) 'end'
            write(13,*) 'title(',"'",'Deformed',"'",')'
            write(13,*) 'axis equal;'
            write(13,*) 'xaxes = get(gca,',"'",'xlim',"'",');'
            write(13,*) 'yaxes = get(gca,',"'",'ylim',"'",');'
            write(13,*) 'axis off;'
            write(13,*) 'subplot(2,1,1)'
            write(13,*) 'hold on'
            write(13,*) 'for e = 1:size(IX,1)'
            write(13,*) '    xx = X(IX(e,1:4),1);'
            write(13,*) '    yy = X(IX(e,1:4),2);'
            write(13,*) '    patch(xx,yy,[1 1 0]);'
            write(13,*) 'end'
            write(13,*) 'title(',"'",'Undeformed',"'",')'
            write(13,*) 'axis([min(xaxes(1),min(X(:,1))) max(xaxes(2),max(X(:,1)))...'
            write(13,*) ' min(yaxes(1),min(X(:,2))) max(yaxes(2),max(X(:,2))) ]);'
            write(13,*) 'axis equal;'
            write(13,*) 'axis off;'
        else
            write (*,'("Unsupported element type in Matlab routine:")')
            write (*,'("plotmatlab -> plotmatlabdef")')
        end if

        ! End file and close
        write(13,*) 'hold off'
        write(13,*) 'set(gcf,',"'",'color',"'",',[ 1  1 1]);'
        close(13)

        ! Print to screen:
        print*,'Matlab plotting called: un/deformed'
        print*,'Generated data file: ', trim(fdata)
        print*,'and script file: ', trim(fscript)
        print*,' '
    end subroutine plotmatlabdef
!
!--------------------------------------------------------------------------------------------------
!
    subroutine plotmatlabeval(title,data1)

        ! Subroutine to plot the element values using Matlab

        character(len=*), intent(in) :: title
        real(wp), dimension(:), intent(in) :: data1
        character(len=100) :: fdata, fscript
        integer i, j, e

        ! write datafile
        imatlab = imatlab + 1
        write (fdata, '(a,"_plotelements_data_",i0.4,".m")') trim(filename), imatlab
        open (13, file=trim(fdata))

        ! write nodal coordinates
        write(13,'("X = [")')
        do i = 1, size(x,1)
            write (13,'(3(f15.9,1x))') (x(i,j), j=1,size(x,2))
        end do
        write(13,'("];")')
        write(13,'( )')

        ! write topology matrix
        write(13,'("IX = [")')
        do e = 1, ne
            write (13,'(20(i6,1x))') (element(e)%ix(i) ,i=1,size(element(e)%ix,1))
        end do
        write(13,'("];")')

        ! write data1/element values
        write(13,'("plotval = [")')
        do i = 1, ne
            write (13,'(f32.15)') data1(i)
        end do
        write(13,'("];")')
        close(13)

        ! Create matlab script
        write(fscript, '(a,"_plotelements_",i0.4,".m")' ) trim(filename), imatlab
        open (13, file = trim(fscript))
        write(13,*) '% Plotting Element Values'
        write(13,*) 'close all'
        write(13,*) 'clear all'
        write(13,*) fdata(1:len_trim(fdata)-2) // ';'
        write(13,*) '% Determine colorscale'
        write(13,*) 'colormap(',"'",'jet',"'",')'
        write(13,*) 'cmap = jet;'
        write(13,*) 'cinterp = linspace(min(plotval),max(plotval),size(cmap,1));'
        write(13,*) '% Make plot'
        write(13,*) 'title(',"'", trim(title) ,"'",')'
        write(13,*) 'hold on'
        ! Element dependent code
        ! NOTE: not possible to mix element types !!!
        if (element(1)%id == 1) then
            write(13,*) 'for e = 1:size(IX,1)'
            write(13,*) '   [dummy,arr_pos] = min(abs(cinterp-plotval(e)));'
            write(13,*) '   xx = X(IX(e,1:2),1);'
            write(13,*) '   yy = X(IX(e,1:2),2);'
            write(13,*) '   plot(xx,yy,',"'",'Color',"'",',cmap(arr_pos,:),',"'",'Linewidth',"'",',1.5);'
            write(13,*) 'end'
        else if (element(1)%id == 2) then
            write(13,*) 'for e = 1:size(IX,1)'
            write(13,*) '   [dummy,arr_pos] = min(abs(cinterp-plotval(e)));'
            write(13,*) '   xx = X(IX(e,1:4),1);'
            write(13,*) '   yy = X(IX(e,1:4),2);'
            write(13,*) '   patch(xx,yy,cmap(arr_pos,:));'
            write(13,*) 'end'
        else
            write (*,'("Unsupported element type in Matlab routine:")')
            write (*,'("plot -> plotmatlabelements")')
        end if

        ! End file and close
        write(13,*) 'axis equal;'
        write(13,*) 'axis off;'
        write(13,*) 'caxis([min(plotval) max(plotval)]);'
        write(13,*) 'colorbar;'
        write(13,*) 'hold off'
        write(13,*) 'set(gcf,',"'",'color',"'",',[ 1  1 1]);'
        close(13)

        ! Print to screen:
        print*,'Matlab plotting called: elements'
        print*,'Generated data file: ', trim(fdata)
        print*,'and script file: ', trim(fscript)
        print*,' '
    end subroutine plotmatlabeval
!
!--------------------------------------------------------------------------------------------------
!
    subroutine plotmatlabevec(title,ppp1,ppp2,pppang)

        ! Vector field plot using Matlab

        real(wp), dimension(:), intent(in) :: ppp1, ppp2, pppang
        character(len=*), intent(in) :: title
        character(len=100) :: fdata, fscript
        integer i, j, e

        ! write datafile
        imatlab = imatlab + 1
        write (fdata, '(a,"_plotvector_data",i0.4,".m")') trim(filename), imatlab
        open (13, file=trim(fdata))

        ! write nodal coordinates
        write(13,'("X = [")')
        do i = 1, size(x,1)
            write (13,'(3(f15.9,1x))') (x(i,j), j=1,size(x,2))
        end do
        write(13,'("];")')
        write(13,'( )')

        ! write topology matrix
        write(13,'("IX = [")')
        do e = 1, ne
            write (13,'(20(i6,1x))') (element(e)%ix(i) ,i=1,size(element(e)%ix,1))
        end do
        write(13,'("];")')

        ! write data1/element values
        write(13,'("vect = [")')
        do i = 1, ne
            write(13,*) ppp1(i), ppp2(i), pppang(i)
        end do
        write(13,'("];")')
        close(13)

        ! Create matlab script
        write(fscript, '(a,"_plotvector_",i0.4,".m")' ) trim(filename), imatlab
        open (13, file = trim(fscript))
        write(13,*) '% Plotting Vector Field, i.e. principle stresses'
        write(13,*) 'close all'
        write(13,*) 'clear all'
        write(13,*) fdata(1:len_trim(fdata)-2) // ';'
        write(13,*) '% Make plot'
        ! Element dependent code
        ! NOTE: not possible to mix element types !!!
        if (element(1)%id == 1) then
            print*,'LINK1 error: cannot do vector plot of truss structure'
            print*,'Files created are empty !!'
        else if (element(1)%id == 2) then
            write(13,*) '% User scale parameter: See fedata - scale_vec'
            write(13,*) 'scale_vec = 1;'
            write(13,*) '% Define characteristic length'
            write(13,*) 'clx = 0;   cly = 0;'
            write(13,*) 'for e = 1:size(IX,1)'
            write(13,*) '    for i = 1:4'
            write(13,*) '        for j = i+1:4'
            write(13,*) '            clx = max(abs(  X(IX(e,i),1) - X(IX(e,j),1)   ),clx);'
            write(13,*) '            cly = max(abs(  X(IX(e,i),2) - X(IX(e,j),2)   ),cly);'
            write(13,*) '        end'
            write(13,*) '    end'
            write(13,*) 'end'
            write(13,*) 'clmax = max(clx, cly);'
            write(13,*) 'scal = max(max(abs(vect(:,1:2))))*sqrt(10)/clmax / scale_vec;'
            write(13,*) '% Make plot'
            write(13,*) 'figure'
            write(13,*) 'hold on'
            write(13,*) 'for e = 1:size(IX,1)'
            write(13,*) '    xx = X(IX(e,1:4),1);'
            write(13,*) '    yy = X(IX(e,1:4),2);'
            write(13,*) '    patch(xx,yy,[1 1 1]);'
            write(13,*) '    % Find approx. center for isoparametric element'
            write(13,*) '    xc = sum(X(IX(e,:),1))/size(IX,2);'
            write(13,*) '    yc = sum(X(IX(e,:),2))/size(IX,2);'
            write(13,*) '    % Directions for vect(:)'
            write(13,*) '    vec = [cos(-vect(e,3)) sin(-vect(e,3)) ...'
            write(13,*) '        cos(-vect(e,3)+pi/2) sin(-vect(e,3)+pi/2)];'
            write(13,*) '    % Plot magnitude and direction of vect_1'
            write(13,*) '    cc = ',"'", 'b',"'",';'
            write(13,*) '    if vect(e,1) < 0,    cc = ',"'",'r',"'",';     end'
            write(13,*) '    quiver(xc,yc,vec(1),vec(2),abs(vect(e,1))/scal,cc)'
            write(13,*) '    quiver(xc,yc,-vec(1),-vec(2),abs(vect(e,1))/scal,cc)'
            write(13,*) '    % Plot magnitude and direction of vect_2'
            write(13,*) '    cc = ',"'",'b',"'",';'
            write(13,*) '    if vect(e,2) < 0,    cc = ',"'",'r',"'",';     end'
            write(13,*) '    quiver(xc,yc,vec(3),vec(4),abs(vect(e,2))/scal,cc)'
            write(13,*) '    quiver(xc,yc,-vec(3),-vec(4),abs(vect(e,2))/scal,cc)'
            write(13,*) 'end   '
        else
            write (*,'("Unsupported element type in Matlab routine:")')
            write (*,'("plot -> plotmatlabevec")')
        endif

        ! End file and close
        write(13,*) 'title( ',"'", trim(title),"'",')'
        write(13,*) 'axis equal;  axis off;  hold off'
        write(13,*) 'set(gcf,',"'",'color',"'",',[ 1  1  1]);'
        close(13)

        ! Print to screen:
        print*,'Matlab plotting called: Vector'
        print*,'Generated data file: ', trim(fdata)
        print*,'and script file: ', trim(fscript)
        print*,' '
    end subroutine plotmatlabevec
!
!--------------------------------------------------------------------------------------------------
!
    subroutine plotmatlabeig(title,freq,evec,times)

        ! Eigenmode plot using Matlab
        ! freq  = eigenvalue
        ! evec  = eigenvector
        ! times = [totaltime, timeinterval]

        real(wp), dimension(:), intent(in) :: evec,times
        real(wp), intent(in) :: freq
        character(len=*), intent(in) :: title
        character(len=100) :: fdata, fscript
        integer i, j, e

        ! write datafile
        imatlab = imatlab + 1
        write (fdata, '(a,"_ploteig_data_",i0.4,".m")') trim(filename), imatlab
        open (13, file=trim(fdata))

        ! write nodal coordinates
        write(13,'("X = [")')
        do i = 1, size(x,1)
            write (13,'(3(f15.9,1x))') (x(i,j), j=1,size(x,2))
        end do
        write(13,'("];")')
        write(13,'( )')

        ! write topology matrix
        write(13,'("IX = [")')
        do e = 1, ne
            write (13,'(20(i6,1x))') (element(e)%ix(i) ,i=1,size(element(e)%ix,1))
        end do
        write(13,'("];")')

        ! write deformation vector
        write(13,'("D = [")')
        do i = 1, neqn
            write (13,'(f15.9,1x)') (evec(i) )
        end do
        write(13,'("];")')
        close(13)

        ! Plot script
        ! Create matlab script
        write(fscript, '(a,"_ploteig_",i0.4,".m")' ) trim(filename), imatlab
        open (13, file = trim(fscript))
        write(13,*) '% Plotting Eigenmodes'
        write(13,*) 'close all'
        write(13,*) 'clear all'
        write(13,*) fdata(1:len_trim(fdata)-2) // ';'
        write(13,*) 'freq = ',freq,';'
        write(13,*) 'timeint = ',times(2),';'
        write(13,*) 'timetot = ',times(1),';'
        write(13,*) '% Find max window size.'
        write(13,*) 'lxmin = min(X(:,1));        lxmax = max(X(:,1));'
        write(13,*) 'lymin = min(X(:,2));        lymax = max(X(:,2));'
        write(13,*) 'dxmin = min(D(1:2:end));    dxmax = max(D(1:2:end));'
        write(13,*) 'dymin = min(D(2:2:end));    dymax = max(D(2:2:end));'
        write(13,*) 'lxmin = lxmin - max(abs(dxmin),abs(dxmax))*1.05;'
        write(13,*) 'lxmax = lxmax + max(abs(dxmin),abs(dxmax))*1.05;'
        write(13,*) 'lymin = lymin - max(abs(dymin),abs(dymax))*1.05;'
        write(13,*) 'lymax = lymax + max(abs(dymin),abs(dymax))*1.05;'
        write(13,*) '% Make plot'
        write(13,*) 'figure'
        write(13,*) 'set(gcf,',"'",'color',"'",',[ 1  1 1]);'
        write(13,*) 'times = 0:timeint:timetot;'
        write(13,*) 'for i = 1:length(times)'
        write(13,*) '    tfact = sin(freq*times(i));'
        write(13,*) '    clf;'
        write(13,*) '    hold on'
        write(13,*) '    for e = 1:size(IX,1)'
        ! Element dependent code
        ! NOTE: not possible to mix element types !!!
        if (element(1)%id == 1) then
            write(13,*) '       edof = [2*IX(e,1)-1 2*IX(e,1) 2*IX(e,2)-1 2*IX(e,2)];'
            write(13,*) '       xx = X(IX(e,1:2),1) + tfact*D(edof(1:2:4));'
            write(13,*) '       yy = X(IX(e,1:2),2) + tfact*D(edof(2:2:4));'
            write(13,*) '       plot(xx,yy,',"'",'b',"',","'",'LineWidth',"'",',1.5)'
        elseif (element(1)%id == 2) then
            write(13,*) '       edof = [2*IX(e,1)-1 2*IX(e,1) 2*IX(e,2)-1 2*IX(e,2)...'
            write(13,*) '          2*IX(e,3)-1 2*IX(e,3) 2*IX(e,4)-1 2*IX(e,4)];'
            write(13,*) '       xx = X(IX(e,1:4),1) + tfact*D(edof(1:2:8));'
            write(13,*) '       yy = X(IX(e,1:4),2) + tfact*D(edof(2:2:8));'
            write(13,*) '       patch(xx,yy,[1 1 0]);'
        else
            write (*,'("Unsupported element type in Matlab routine:")')
            write (*,'("plotmatlab -> plotmatlabdef")')
        endif

        ! End file and close
        write(13,*) '    end'
        write(13,*) '    axis([lxmin lxmax lymin lymax])'
        write(13,*) '    axis off'
        write(13,*) '    title( ',"'", trim(title),"'",')'
        write(13,*) '    pause(0.01)'
        write(13,*) 'end'
        close(13)

        ! Print to screen:
        print*,'Matlab plotting called: eigenmode'
        print*,'Generated data file: ', trim(fdata)
        print*,'and script file: ', trim(fscript)
        print*,' '
    end subroutine plotmatlabeig

end module processor
