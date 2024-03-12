module search
    use interactions
    use domain
    use functions
    use initialize
    use kernel


    implicit none
    
    contains
    
    subroutine nsearch2(list,n,celln,cell1,cell2,cell3,cell4,smln)
        implicit none

        type(verlet),intent(inout) :: list
        type(cell),intent(inout),target :: celln,cell1,cell2,cell3,cell4
        real(dp),intent(in) :: smln
        integer,intent(in) :: n
        integer ::counter2
        real(dp) :: radial
        real(dp) ::res2,res1,res3,res4,res5,res6,t1,t2
        real(dp),dimension(4) :: res,res_s
        real(dp) :: heff
        integer :: tempcal


        list%count=0

        res=0.0_dp
        res2=0.0_dp
        radial=0.0_dp
        res1=0.0_dp
        res3=0.0_dp
        res4=0.0_dp

        res5=0.0_dp
        res6=0.0_dp
        res_s=0.0_dp

        celln%plist(n)%free=.false.
        celln%pplist(n)%gradvx=0.0_dp
        t1=0.0_dp
        t2=0.0_dp

        tempcal=0

        heff=smln

        do counter2=1,cell1%ptot
            
            
                    if ((celln%plist(n)%tid/=cell1%plist(counter2)%tid) &
                    .or.((celln%plist(n)%tid==cell1%plist(counter2)%tid) &
                    .and.(celln%plist(n)%matid/=cell1%plist(counter2)%matid))) then

                    radial=dist(celln%plist(n),cell1%plist(counter2))

                    if((radial<(2*smln)))then

                    list%count=list%count+1
                    ! list%interlist(1,list%count)=cell1%cellid(1)
                    ! list%interlist(2,list%count)=cell1%cellid(2)
                    ! list%interlist(3,list%count)=counter2
                    list%nh(list%count)%part=>cell1%plist(counter2)
                    list%pnh(list%count)%ppart=>cell1%pplist(counter2)
                    list%dist(list%count)=radial

                    !Gradient Correction

                    res2=Wabx(cell1%plist(counter2),celln%plist(n),radial,smln)  
                    res4=Waby(cell1%plist(counter2),celln%plist(n),radial,smln)
                    list%klt(1,list%count)=res2
                    list%klt(2,list%count)=res4 
                    ! list%klt(3,list%count)=Wo(smln)    

                    res1=(cell1%plist(counter2)%mass)*(cell1%plist(counter2)%x-celln%plist(n)%x) &
                    /(cell1%plist(counter2)%density)
                    res3=(cell1%plist(counter2)%mass)*(cell1%plist(counter2)%y-celln%plist(n)%y) &
                    /(cell1%plist(counter2)%density)

                    res(1)=res(1)+res2*res1
                    res(2)=res(2)+res4*res1
                    res(3)=res(3)+res2*res3
                    res(4)=res(4)+res4*res3

                    ! res5=res5+(cell1%plist(counter2)%mass)*(cell1%plist(counter2)%x-celln%plist(n)%x) &
                    ! *Wab(radial,smln)/(cell1%plist(counter2)%density)
                    ! res6=res6+(cell1%plist(counter2)%mass)*(cell1%plist(counter2)%y-celln%plist(n)%y) &
                    ! *Wab(radial,smln)/(cell1%plist(counter2)%density)

                    ! celln%pplist(n)%shep=celln%pplist(n)%shep+ &
                    ! (Wab(radial,smln)*(cell1%plist(counter2)%mass)/(cell1%plist(counter2)%density))

                    ! celln%pplist(n)%shep_gradx=celln%pplist(n)%shep_gradx+ &
                    ! (res2*(cell1%plist(counter2)%mass)/(cell1%plist(counter2)%density))

                    ! celln%pplist(n)%shep_grady=celln%pplist(n)%shep_grady+ &
                    ! (res4*(cell1%plist(counter2)%mass)/(cell1%plist(counter2)%density))

                    end if

                end if
        end do

        do counter2=1,cell2%ptot
            
                    if ((celln%plist(n)%tid/=cell2%plist(counter2)%tid) &
                    .or.((celln%plist(n)%tid==cell2%plist(counter2)%tid) &
                    .and.(celln%plist(n)%matid/=cell2%plist(counter2)%matid))) then

                    radial=dist(celln%plist(n),cell2%plist(counter2))
                    if((radial<(2*smln)))then

                    list%count=list%count+1
                    ! list%interlist(1,list%count)=cell2%cellid(1)
                    ! list%interlist(2,list%count)=cell2%cellid(2)
                    ! list%interlist(3,list%count)=counter2
                    list%nh(list%count)%part=>cell2%plist(counter2)
                    list%pnh(list%count)%ppart=>cell2%pplist(counter2)
                    list%dist(list%count)=radial

                    !Gradient Correction

                    res2=Wabx(cell2%plist(counter2),celln%plist(n),radial,smln)  
                    res4=Waby(cell2%plist(counter2),celln%plist(n),radial,smln)
                    list%klt(1,list%count)=res2
                    list%klt(2,list%count)=res4 
                    ! list%klt(3,list%count)=Wo(smln)    

                    res1=(cell2%plist(counter2)%mass)*(cell2%plist(counter2)%x-celln%plist(n)%x) &
                    /(cell2%plist(counter2)%density)
                    res3=(cell2%plist(counter2)%mass)*(cell2%plist(counter2)%y-celln%plist(n)%y) &
                    /(cell2%plist(counter2)%density)

                    res(1)=res(1)+res2*res1
                    res(2)=res(2)+res4*res1
                    res(3)=res(3)+res2*res3
                    res(4)=res(4)+res4*res3

                    ! res5=res5+(cell2%plist(counter2)%mass)*(cell2%plist(counter2)%x-celln%plist(n)%x) &
                    ! *Wab(radial,smln)/(cell2%plist(counter2)%density)
                    ! res6=res6+(cell2%plist(counter2)%mass)*(cell2%plist(counter2)%y-celln%plist(n)%y) &
                    ! *Wab(radial,smln)/(cell2%plist(counter2)%density)

                    ! celln%pplist(n)%shep=celln%pplist(n)%shep+ &
                    ! (Wab(radial,smln)*(cell2%plist(counter2)%mass)/(cell2%plist(counter2)%density))

                    ! celln%pplist(n)%shep_gradx=celln%pplist(n)%shep_gradx+ &
                    ! (res2*(cell2%plist(counter2)%mass)/(cell2%plist(counter2)%density))

                    ! celln%pplist(n)%shep_grady=celln%pplist(n)%shep_grady+ &
                    ! (res4*(cell2%plist(counter2)%mass)/(cell2%plist(counter2)%density))

                    end if

            end if
        end do

        do counter2=1,cell3%ptot
            
                    if ((celln%plist(n)%tid/=cell3%plist(counter2)%tid) &
                    .or.((celln%plist(n)%tid==cell3%plist(counter2)%tid) &
                    .and.(celln%plist(n)%matid/=cell3%plist(counter2)%matid))) then

                    radial=dist(celln%plist(n),cell3%plist(counter2))
                    if((radial<(2*smln)))then

                    list%count=list%count+1
                    ! list%interlist(1,list%count)=cell3%cellid(1)
                    ! list%interlist(2,list%count)=cell3%cellid(2)
                    ! list%interlist(3,list%count)=counter2
                    list%nh(list%count)%part=>cell3%plist(counter2)
                    list%pnh(list%count)%ppart=>cell3%pplist(counter2)
                    list%dist(list%count)=radial

                    !Gradient Correction

                    res2=Wabx(cell3%plist(counter2),celln%plist(n),radial,smln)  
                    res4=Waby(cell3%plist(counter2),celln%plist(n),radial,smln)
                    list%klt(1,list%count)=res2
                    list%klt(2,list%count)=res4 
                    ! list%klt(3,list%count)=Wo(smln)    

                    res1=(cell3%plist(counter2)%mass)*(cell3%plist(counter2)%x-celln%plist(n)%x) &
                    /(cell3%plist(counter2)%density)
                    res3=(cell3%plist(counter2)%mass)*(cell3%plist(counter2)%y-celln%plist(n)%y) &
                    /(cell3%plist(counter2)%density)

                    res(1)=res(1)+res2*res1
                    res(2)=res(2)+res4*res1
                    res(3)=res(3)+res2*res3
                    res(4)=res(4)+res4*res3

                    ! res5=res5+(cell3%plist(counter2)%mass)*(cell3%plist(counter2)%x-celln%plist(n)%x) &
                    ! *Wab(radial,smln)/(cell3%plist(counter2)%density)
                    ! res6=res6+(cell3%plist(counter2)%mass)*(cell3%plist(counter2)%y-celln%plist(n)%y) &
                    ! *Wab(radial,smln)/(cell3%plist(counter2)%density)

                    ! celln%pplist(n)%shep=celln%pplist(n)%shep+ &
                    ! (Wab(radial,smln)*(cell3%plist(counter2)%mass)/(cell3%plist(counter2)%density))

                    ! celln%pplist(n)%shep_gradx=celln%pplist(n)%shep_gradx+ &
                    ! (res2*(cell3%plist(counter2)%mass)/(cell3%plist(counter2)%density))

                    ! celln%pplist(n)%shep_grady=celln%pplist(n)%shep_grady+ &
                    ! (res4*(cell3%plist(counter2)%mass)/(cell3%plist(counter2)%density))

                    end if

            end if
        end do

        do counter2=1,cell4%ptot
            
                    if ((celln%plist(n)%tid/=cell4%plist(counter2)%tid) &
                    .or.((celln%plist(n)%tid==cell4%plist(counter2)%tid) &
                    .and.(celln%plist(n)%matid/=cell4%plist(counter2)%matid))) then

                    radial=dist(celln%plist(n),cell4%plist(counter2))
                    if((radial<(2*smln)))then

                    list%count=list%count+1
                    ! list%interlist(1,list%count)=cell4%cellid(1)
                    ! list%interlist(2,list%count)=cell4%cellid(2)
                    ! list%interlist(3,list%count)=counter2
                    list%nh(list%count)%part=>cell4%plist(counter2)
                    list%pnh(list%count)%ppart=>cell4%pplist(counter2)
                    list%dist(list%count)=radial

                    !Gradient Correction

                    res2=Wabx(cell4%plist(counter2),celln%plist(n),radial,smln)  
                    res4=Waby(cell4%plist(counter2),celln%plist(n),radial,smln)
                    list%klt(1,list%count)=res2
                    list%klt(2,list%count)=res4 
                    ! list%klt(3,list%count)=Wo(smln)   

                    res1=(cell4%plist(counter2)%mass)*(cell4%plist(counter2)%x-celln%plist(n)%x) &
                    /(cell4%plist(counter2)%density)
                    res3=(cell4%plist(counter2)%mass)*(cell4%plist(counter2)%y-celln%plist(n)%y) &
                    /(cell4%plist(counter2)%density)

                    res(1)=res(1)+res2*res1
                    res(2)=res(2)+res4*res1
                    res(3)=res(3)+res2*res3
                    res(4)=res(4)+res4*res3

                    ! res5=res5+(cell4%plist(counter2)%mass)*(cell4%plist(counter2)%x-celln%plist(n)%x) &
                    ! *Wab(radial,smln)/(cell4%plist(counter2)%density)
                    ! res6=res6+(cell4%plist(counter2)%mass)*(cell4%plist(counter2)%y-celln%plist(n)%y) &
                    ! *Wab(radial,smln)/(cell4%plist(counter2)%density)

                    ! celln%pplist(n)%shep=celln%pplist(n)%shep+ &
                    ! (Wab(radial,smln)*(cell4%plist(counter2)%mass)/(cell4%plist(counter2)%density))

                    ! celln%pplist(n)%shep_gradx=celln%pplist(n)%shep_gradx+ &
                    ! (res2*(cell4%plist(counter2)%mass)/(cell4%plist(counter2)%density))

                    ! celln%pplist(n)%shep_grady=celln%pplist(n)%shep_grady+ &
                    ! (res4*(cell4%plist(counter2)%mass)/(cell4%plist(counter2)%density))

                    end if

            end if
        end do

        if (celln%plist(n)%tid==3) then
        
        celln%plist(n)%free=.false.
        celln%plist(n)%vicinity=.false.
        celln%pplist(n)%gradvx=(res(1)+res(4))!t1+t2

        res_s=res

        call invertmat2D(res)

        celln%pplist(n)%coff=res

        ! res_s(1)=(res_s(1)/celln%pplist(n)%shep)- &
        !             (res5*celln%pplist(n)%shep_gradx/ &
        !             (celln%pplist(n)%shep**2))

        ! res_s(2)=(res_s(2)/celln%pplist(n)%shep)- &
        !         (res5*celln%pplist(n)%shep_grady/ &
        !         (celln%pplist(n)%shep**2)) 
                
                
        ! res_s(3)=(res_s(3)/celln%pplist(n)%shep)- &
        !         (res6*celln%pplist(n)%shep_gradx/ &
        !         (celln%pplist(n)%shep**2))

        ! res_s(4)=(res_s(4)/celln%pplist(n)%shep)- &
        !         (res6*celln%pplist(n)%shep_grady/ &
        !         (celln%pplist(n)%shep**2)) 

        ! call invertmat2D(res_s)

        ! celln%pplist(n)%gradcoff=res_s

        
        if((celln%pplist(n)%gradvx<=(lamfs*maxdivr)).and.(.not.(celln%plist(n)%buffer))) then
            celln%plist(n)%free=.true.
        end if

        end if

    end subroutine nsearch2

    subroutine neighbour()

        implicit none

        type(cell),pointer :: celln,cell1,cell2,cell3,cell4

        integer :: i,j,k,m

        !$omp do private(k,i,j,cell1,cell2,cell3,cell4,celln) schedule(runtime) collapse(2)
        do j=sx,ex
            do i=sy,ey
            ! if (dpcell(i,j)%ptot/=0) then
    
                do k=1,dpcell(i,j)%ptot

                    celln=>dpcell(i,j)


                if(dpcell(i,j)%plist(k)%x<=(dpcell(i,j)%xleft+ &
                abs((dpcell(i,j)%xleft-dpcell(i,j)%xright)/2)).and. &
                dpcell(i,j)%plist(k)%y<=(dpcell(i,j)%ybot+ &
                abs((dpcell(i,j)%ytop-dpcell(i,j)%ybot)/2))) then
                    ! call nsearch(dpcell(i,j)%list(k),k,dpcell(i,j),dpcell(i,j-1), &
                    ! dpcell(i+1,j),dpcell(i+1,j-1),h1)

                    cell1=>dpcell(i,j-1)
                    cell2=>dpcell(i+1,j-1)
                    cell3=>dpcell(i,j)
                    cell4=>dpcell(i+1,j)

                else if(dpcell(i,j)%plist(k)%x>(dpcell(i,j)%xleft+ &
                abs((dpcell(i,j)%xleft-dpcell(i,j)%xright)/2)).and. &
                dpcell(i,j)%plist(k)%y<=(dpcell(i,j)%ybot+ &
                abs((dpcell(i,j)%ytop-dpcell(i,j)%ybot)/2))) then
                    ! call nsearch(dpcell(i,j)%list(k),k,dpcell(i,j),dpcell(i,j+1), &
                    ! dpcell(i+1,j),dpcell(i+1,j+1),h1)

                    cell1=>dpcell(i,j)
                    cell2=>dpcell(i+1,j)
                    cell3=>dpcell(i,j+1)
                    cell4=>dpcell(i+1,j+1)

                else if(dpcell(i,j)%plist(k)%x>(dpcell(i,j)%xleft+ &
                abs((dpcell(i,j)%xleft-dpcell(i,j)%xright)/2)).and. &
                dpcell(i,j)%plist(k)%y>(dpcell(i,j)%ybot+ &
                abs((dpcell(i,j)%ytop-dpcell(i,j)%ybot)/2))) then
                    ! call nsearch(dpcell(i,j)%list(k),k,dpcell(i,j),dpcell(i-1,j+1), &
                    ! dpcell(i-1,j),dpcell(i,j+1),h1)

                    cell1=>dpcell(i-1,j)
                    cell2=>dpcell(i,j)
                    cell3=>dpcell(i-1,j+1)
                    cell4=>dpcell(i,j+1)

                else if(dpcell(i,j)%plist(k)%x<=(dpcell(i,j)%xleft+ &
                abs((dpcell(i,j)%xleft-dpcell(i,j)%xright)/2)).and. &
                dpcell(i,j)%plist(k)%y>(dpcell(i,j)%ybot+ &
                abs((dpcell(i,j)%ytop-dpcell(i,j)%ybot)/2))) then
                    ! call nsearch(dpcell(i,j)%list(k),k,dpcell(i,j),dpcell(i-1,j), &
                    ! dpcell(i-1,j-1),dpcell(i,j-1),h1)

                    cell1=>dpcell(i-1,j-1)
                    cell2=>dpcell(i,j-1)
                    cell3=>dpcell(i-1,j)
                    cell4=>dpcell(i,j)

                else if(dpcell(i,j)%plist(k)%x==(dpcell(i,j)%xleft+ &
                    abs((dpcell(i,j)%xleft-dpcell(i,j)%xright)/2)).and. &
                    dpcell(i,j)%plist(k)%y==(dpcell(i,j)%ybot+ &
                    abs((dpcell(i,j)%ytop-dpcell(i,j)%ybot)/2))) then
                    ! call nsearch(dpcell(i,j)%list(k),k,dpcell(i,j),dpcell(i-1,j), &
                    ! dpcell(i-1,j-1),dpcell(i,j-1),h1)

                    cell1=>dpcell(i-1,j-1)
                    cell2=>dpcell(i,j-1)
                    cell3=>dpcell(i-1,j)
                    cell4=>dpcell(i,j)

                end if

                call nsearch2(dpcell(i,j)%list(k),k,celln,cell1,cell2,cell3,cell4,h1)

                end do
            ! end if
            end do               
        end do
        !$omp end do

        return
    end subroutine neighbour

    
end module search