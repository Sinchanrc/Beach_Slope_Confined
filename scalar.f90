module scalar

    use particle
    use initialize
    use domain
    use kernel

    implicit none


    contains

    subroutine dispersivity(dif,den,poro,vel,al,at)
        real(dp), intent(in) :: dif,den,poro,vel
        real(dp), intent(out) ::  al,at
        real(dp) :: Pec,Sch,tor,Dmt,Re,D1l,D1t
        integer :: case
    
        tor=1.106_dp*poro**(-0.196_dp)
        Pec=vel*0.00102_dp/(dif*poro)
        Pec=Pec*tor
        Sch=mu/(den*dif)
        Dmt=dif/tor
        Re=den*vel*0.00102_dp/(mu*poro)
    
        if (Pec<0.1_dp) then
        case=1
        elseif ((Pec>=0.1_dp).and.(Pec<4.0_dp)) then
        case=2
        elseif((Pec>=4.0_dp).and.(Re<10.0_dp)) then
        case=3
        elseif ((Re>=10.0_dp).and.(Pec<1e6)) then
        case=4
        end if
    
        select case (case)
        case(1)
            D1l=Dmt
        case(2)
            D1l=(Pec/(0.8_dp*Pec+0.4_dp))*Dmt
        case(3)
            D1l=(Pec/sqrt((18.0_dp*Pec**(-1.2_dp))+2.35_dp*Sch**(-0.38_dp)))*Dmt
        case(4)
            D1l=(Pec/(((25.0_dp*Sch**(1.14_dp))/Pec)+0.5_dp))*Dmt
        end select
    
        al=D1l*poro/vel
    
        if (Pec<1.0_dp) then
            case=1
        elseif ((Pec>=1.0_dp).and.(Pec<1600.0_dp)) then
            if (Sch<550.0_dp) then
                case=2
                else
                case=3
            end if
        elseif((Pec>=1600.0_dp).and.(Pec<1e6)) then
            if (Sch<550.0_dp) then
                case=4
                else
                case=5
            end if
        end if
    
        select case (case)
        case(1)
            D1t=Dmt
        case(2)
            D1t=(1.0_dp+1.0_dp/((2.7_dp*1e-5*Sch)+(12.0_dp/Pec)))*Dmt
        case(3)
            D1t=(1.0_dp+1.0_dp/(0.017_dp+(14.0_dp/Pec)))*Dmt
        case(4)
            D1t=(Pec/((0.058_dp*Sch+14.0_dp)-((0.058_dp*Sch+2.0_dp) &
            *exp((-500.0_dp*Sch**(0.5_dp))/Pec))))*Dmt
        case(5)
            D1t=(Pec/(45.9_dp-(33.9_dp*exp(-21.0_dp*Sch/Pec))))*Dmt
    
        end select
    
        at=D1t*poro/vel
    
        
    end subroutine

    subroutine massupdate(stime)

        use particle
    
        implicit none 
        integer :: i,j,k
        real(dp),intent(in) :: stime
    
        !$omp do schedule (runtime) private(k,i,j) collapse(2)
            do j=sx,ex
                do i=sy,ey            
                    do k=1,dpcell(i,j)%ptot
    
                    if((dpcell(i,j)%plist(k)%tid==3)) then
    
    
                        dpcell(i,j)%plist(k)%con=dpcell(i,j)%plist(k)%con+ &
                        dpcell(i,j)%pplist(k)%cdiff*stime
    
    
                    end if
    
                    end do
                end do
            end do
        !$omp end do
    
    end subroutine

    subroutine scalart

        integer :: i,j,k,m
        real(dp) :: t1x,t1y,t2,t3,con1,con2,Dxi,Dyi,Dxj,Dyj,vi,vj
        real(dp),parameter :: al=0.004_dp,at=0.0004_dp,Dm=1e-9,tor=1.0_dp
        real(dp) :: ali,ati,alj,atj

        !$omp do schedule (runtime) collapse(2) &
        !$omp private(m,k,i,j,t1x,t1y,t2,t3,con1,con2,Dxi,Dyi,Dxj,Dyj,vi,vj,ali,ati,alj,atj)
            do j=sx,ex
                do i=sy,ey            
            
                do k=1,dpcell(i,j)%ptot

                if((dpcell(i,j)%plist(k)%tid==3)) then

                dpcell(i,j)%pplist(k)%cdiff=0.0_dp

                con1=dpcell(i,j)%plist(k)%con

                vi=sqrt(dpcell(i,j)%plist(k)%vx**2+dpcell(i,j)%plist(k)%vy**2)

                ! call dispersivity(Dm,dpcell(i,j)%plist(k)%oden,&
                ! dpcell(i,j)%pplist(k)%porosity,vi,ali,ati)

                ! Dxi=Dm+((ali*dpcell(i,j)%plist(k)%vx**2)+(ati*dpcell(i,j)%plist(k)%vy**2)+ &
                ! ((ali-ati)*dpcell(i,j)%plist(k)%vx*dpcell(i,j)%plist(k)%vy))/(vi*dpcell(i,j)%pplist(k)%porosity)

                Dxi=Dm+((al*dpcell(i,j)%plist(k)%vx**2)+(at*dpcell(i,j)%plist(k)%vy**2)+ &
                ((al-at)*dpcell(i,j)%plist(k)%vx*dpcell(i,j)%plist(k)%vy))/(vi*dpcell(i,j)%pplist(k)%porosity)

                    Dxi=merge(Dxi,Dm,Dxi>0.0_dp)
                    Dxi=merge(Dxi,0.12_dp*dl**2/dt,Dxi<=0.12_dp*dl**2/dt)

                ! Dyi=Dm+((ati*dpcell(i,j)%plist(k)%vx**2)+(ali*dpcell(i,j)%plist(k)%vy**2)+ &
                ! ((ali-ati)*dpcell(i,j)%plist(k)%vx*dpcell(i,j)%plist(k)%vy))/(vi*dpcell(i,j)%pplist(k)%porosity)

                Dyi=Dm+((at*dpcell(i,j)%plist(k)%vx**2)+(al*dpcell(i,j)%plist(k)%vy**2)+ &
                ((al-at)*dpcell(i,j)%plist(k)%vx*dpcell(i,j)%plist(k)%vy))/(vi*dpcell(i,j)%pplist(k)%porosity)

                Dyi=merge(Dyi,Dm,Dyi>0.0_dp)
                Dyi=merge(Dyi,0.12_dp*dl**2/dt,Dyi<=0.12_dp*dl**2/dt)

                do m=1,dpcell(i,j)%list(k)%count
                    associate(x=>dpcell(i,j)%list(k)%nh(m)%part, &
                        y=>dpcell(i,j)%list(k)%pnh(m)%ppart, &
                        z=>dpcell(i,j)%list(k)%klt)

                    if ((x%tid==3)) then

                        con2=x%con

                        vj=sqrt(x%vx**2+x%vy**2)

                        ! call dispersivity(Dm,x%oden,y%porosity,vj,alj,atj)

                        ! Dxj=Dm+((alj*x%vx**2)+(atj*x%vy**2)+((alj-atj)*x%vx*x%vy))/(vj*y%porosity)

                        ! Dyj=Dm+((atj*x%vx**2)+(alj*x%vy**2)+((alj-atj)*x%vx*x%vy))/(vj*y%porosity)

                        Dxj=Dm+((al*x%vx**2)+(at*x%vy**2)+((al-at)*x%vx*x%vy))/(vj*y%porosity)

                        Dyj=Dm+((at*x%vx**2)+(al*x%vy**2)+((al-at)*x%vx*x%vy))/(vj*y%porosity)

                        Dxj=merge(Dxj,Dm,Dxj>0.0_dp)
                        Dxj=merge(Dxj,0.12_dp*dl**2/dt,Dxj<=0.12_dp*dl**2/dt)

                        Dyj=merge(Dyj,Dm,Dyj>0.0_dp)
                        Dyj=merge(Dyj,0.12_dp*dl**2/dt,Dyj<=0.12_dp*dl**2/dt)

                        if ((dpcell(i,j)%pplist(k)%porosity>=0.5_dp).or.(y%porosity>=0.5_dp)) then

                        if ((dpcell(i,j)%pplist(k)%nut<1e-6).and.(y%nut<1e-6)) then

                            t1x=Dm+(dpcell(i,j)%pplist(k)%nut+y%nut)*0.50_dp/tschmidt
                        else

                        t1x=2*(((Dm+dpcell(i,j)%pplist(k)%nut/tschmidt)*(Dm+y%nut/tschmidt))/&
                        ((Dm+dpcell(i,j)%pplist(k)%nut/tschmidt)+(Dm+y%nut/tschmidt)))
                        

                        end if
                        t1y=t1x

                        else

                            t1x=2.0_dp*tor*Dxi*Dxj/(Dxi+Dxj)
                            t1y=2.0_dp*tor*Dyi*Dyj/(Dyi+Dyj)


                        end if

                        ! t1=Dm

                        t2=2*dpcell(i,j)%pplist(k)%porosity*y%porosity/&
                            (dpcell(i,j)%pplist(k)%porosity+y%porosity)
                        
                        t3=(dpcell(i,j)%plist(k)%ovol)*(dpcell(i,j)%pplist(k)%porosity+y%porosity)/&
                        (dpcell(i,j)%pplist(k)%porosity*y%porosity)

                        ! dpcell(i,j)%pplist(k)%cdiff=dpcell(i,j)%pplist(k)%cdiff + &
                        ! (t1x*t2*t3*(con1-con2)* &
                        ! (dpcell(i,j)%plist(k)%x-x%x)*z(1,m)/(dpcell(i,j)%list(k)%dist(m)**2+lam))+ &
                        ! (t1y*t2*t3*(con1-con2)* &
                        ! (dpcell(i,j)%plist(k)%y-x%y)*z(2,m)/(dpcell(i,j)%list(k)%dist(m)**2+lam))

                        dpcell(i,j)%pplist(k)%cdiff=dpcell(i,j)%pplist(k)%cdiff + &
                        (t1x*t2*t3*(con1-con2)* &
                        (dpcell(i,j)%plist(k)%x-x%x)*hWabx(x,dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),h1)/(dpcell(i,j)%list(k)%dist(m)**2+lam))+ &
                        (t1y*t2*t3*(con1-con2)* &
                        (dpcell(i,j)%plist(k)%y-x%y)*hWaby(x,dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),h1)/(dpcell(i,j)%list(k)%dist(m)**2+lam))

                    end if
                    end associate

                end do

                end if

                end do

                end do
            end do
        !$omp end do

    end subroutine

    subroutine scalarupdate(stime)

        use particle

        implicit none 
        integer :: i,j,k
        real(dp),intent(in) :: stime

        !$omp do schedule (runtime) private(k,i,j) collapse(2)
            do j=sx,ex
                do i=sy,ey            
                    do k=1,dpcell(i,j)%ptot

                    if((dpcell(i,j)%plist(k)%tid==3)) then

                        dpcell(i,j)%plist(k)%con=dpcell(i,j)%plist(k)%con+ &
                        dpcell(i,j)%pplist(k)%cdiff*stime

                        if ((dpcell(i,j)%plist(k)%buffer).and.(dpcell(i,j)%plist(k)%x<1.6_dp)) then

                            dpcell(i,j)%plist(k)%con=0.5_dp

                        end if
                        


                    end if

                    end do
                end do
            end do
        !$omp end do

    end subroutine

    subroutine densityupdate

        implicit none 

        integer :: i,j,k

        !$omp do schedule (runtime) private(k,i,j) collapse(2)
            do j=sx,ex
                do i=sy,ey            
                    do k=1,dpcell(i,j)%ptot

                    if((dpcell(i,j)%plist(k)%tid==3)) then

                        ! dpcell(i,j)%plist(k)%density=0.5*(rhomax+rhomin)+dpcell(i,j)%plist(k)%con*(rhomax-rhomin)
                        ! dpcell(i,j)%plist(k)%mass=dpcell(i,j)%plist(k)%density*dpcell(i,j)%plist(k)%ovol
                        dpcell(i,j)%plist(k)%oden=0.5*(rhomax+rhomin)+dpcell(i,j)%plist(k)%con*(rhomax-rhomin)
                        dpcell(i,j)%plist(k)%mass=dpcell(i,j)%plist(k)%oden*dpcell(i,j)%plist(k)%ovol
                        dpcell(i,j)%plist(k)%density=dpcell(i,j)%plist(k)%oden*dpcell(i,j)%pplist(k)%porosity


                    end if

                    end do
                end do
            end do
        !$omp end do

    end subroutine
    

end module
