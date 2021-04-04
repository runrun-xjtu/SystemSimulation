
	program main

      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      parameter (ncmax=20)   !最大的工质数（应该用不到）
      dimension x(ncmax),xliq(ncmax),xvap(ncmax),f(ncmax)
      character hrf*3, herr*255   !以下三行应该没什么用，不过别删
      character*255 hf(ncmax),hfmix
      dimension xdew(ncmax),ybub(ncmax)
	integer i,j,n
	real h1(100),h2(100),h3(100),h4(100),h5(100),t6(100),h7(100),v1,v3(100),s1(100),s3(100),p1(100),p5(100),p3(100)    !声明各点的焓 温度 比体积 压力等
    real q0(100),wd(100),wg(100),qmd(100),qmg(100),qvd(100),qvg(100)
	real pg(100),pd(100),cop(100),xm(100),xv(100),qk(100),tm(100),pm(100),ha,sa
	real nisd(100),nisg(100),t2(100),t4(100),h6(100),h8(100),s2(100),s4(100),s5(100),s6(100),s7(100),s8(100)
	real x6(100),x8(100),h9(100),s9(100),h10(100),s10(100),ex1(100),ex2(100),ex3(100),ex4(100),ex5(100),ex6(100),ex7(100),ex8(100)
	real exloss_d(100),exloss_g(100),exloss_con(100),exloss_eva(100),exloss_sep(100),exloss_thrg(100),exloss_thrd(100),qeva(100),exloss_back(100)
    real mass_mol,ex_eff(100),dt_backheat,t_hotout,t_coldout,ex77(100),t77(100),h77(100),s77(100),ex1_sat(100),h1_sat(100),s1_sat(100),t1(100),j_dhmin,dh_backheat(101),h_hotout,h_coldout,h_hotin,h_coldin,pmm,pl

    !*****************输入数据******************************   

    real ::tl=243.15, th=308.15, t5=308.15, ta=298.15, pa=101.325,dp=0                       !   T-K   p-kPa
    real ::capacity=2500, ha1_con=300.43, ha2_con=306.46, ha1_eva=247.14, ha2_eva=251.16       !   Q-W   h-kJ/kg 
	real ::ex1_con=0, ex2_con=0.096818, ex1_eva=4.171, ex2_eva=4.951   

	real ::backheat_ratio=2  !回热器冷流体出口与热流体进口温差 K

    hf='r22.fld'
    !*****************输入数据******************************


      nc=1                      !下面5行是程序初始化需要的
      hfmix='hmx.bnc'          
      hrf='DEF'                 
      call SETUP (nc,hf,hfmix,hrf,ierr,herr)

    n=100
	qmd=1                  !以qmd(低压压缩机进气量)为单位“1”
  
	open(1,file='result.txt')
	write(1,10)

    !*****************变工况******************************
	dt=(th-tl-10)/n      
    do i=1,n-1             
    tm(i)=tl+5+dt*i  
	
	kph=2;t=tm(i)     !计算3点中间压力（气态）
	call SATT (t,x,kph,p,rhol,rhov,xliq,xvap,ierr,herr)
	p3(i)=p	
	        
    end do

    !*****************变工况********************************

    !real mass_mol,ex_eff(100),dt_backheat,t_hotout(100),t_coldout(100),t7(100),t11(100),j_dhmin,dh_backheat(100)

    do i=1,n-1

    dh_backheat(1)=100000
	dt_backheat=(tm(i)-tl)/100

    do j=1,100
    t_hotout=tl+j*dt_backheat
	t_coldout=tm(i)-j*dt_backheat/backheat_ratio

	kph=2;t=tl      !中间压力，低压 
    call SATT (t,x,kph,p,rhol,rhov,xliq,xvap,ierr,herr);      
	pl=p;             
	kph=1;t=tm(i)
	call SATT (t,x,kph,p,rhol,rhov,xliq,xvap,ierr,herr)
	pmm=p;

	t=tl+0.1;p=pl    !四点焓
    call TPFLSH (t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
    h_coldin=h
	t=tm(i)-0.1;p=pmm    
    call TPFLSH (t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
    h_hotin=h
	t=t_coldout;p=pl    !四点焓
    call TPFLSH (t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
    h_coldout=h
	t=t_hotout;p=pmm    
    call TPFLSH (t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
    h_hotout=h

    dh_backheat(j+1)=abs((h_hotin-h_hotout)-(h_coldout-h_coldin))
	
	if(dh_backheat(j+1)<dh_backheat(j)) then
    j_dhmin=j

    end if
    end do

	write(*,*),j_dhmin

    t77(i)=tl+j_dhmin*dt_backheat
	t1(i)=tm(i)-j_dhmin*dt_backheat/backheat_ratio

    !write(*,*),'t1*  t1, t7*  t7',t11(i),tl,t7(i),tm(i)
    end do

    !*****************计算回热 准备条件********************************
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
    do i=1,n-1        !状态点计算循环开始   ! ! GENERAL FLASH SUBROUTINES

!    if (i/=50)then    !!!!!!!单工况
!	cycle
!	end if 

	kph=2;t=tl        ! kph = 2代表为气态 kph = 1 代表液态
    call SATT (t,x,kph,p,rhol,rhov,xliq,xvap,ierr,herr);      
	p1(i)=p;             !计算蒸发压力、1点（气态）压力

	kph=1;t=th
	call SATT (t,x,kph,p,rhol,rhov,xliq,xvap,ierr,herr)
	p5(i)=p;             !计算冷凝压力、5点（液态）压力

	t=t1(i);p=p1(i)         !根据1点 p T 计算h  v1 s1
    call TPFLSH (t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
    h1(i)=h;s1(i)=s    ! d--bulk molar density [mol/L] 0.05882 = 1/17 换算后单位为m3/kg         v1=0.05882/d;

	t=t5-0.1;p=p5(i)         !根据5点 p T 计算h5  s5
    call TPFLSH (t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
    h5(i)=h;s5(i)=s 

    t=tl+0.1;p=p1(i)         
    call TPFLSH (t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
    h1_sat(i)=h;s1_sat(i)=s


	! *********以上 1 , 5,1'(11) ****************************

    !p3(i)=sqrt(p1(i)*p5(i))
    !kph=2;p3(i)=p
	!call SATP (p,x,kph,t,rhol,rhov,xliq,xvap,ierr,herr)
	!tm(i)=t
	
	! *********以上 p3 *******************************

    nisd(i)=0.874-0.0135*(p3(i)/p1(i))    !等熵效率
	nisg(i)=0.874-0.0135*(p5(i)/p3(i))    !等熵效率
	write(*,*),'等s效率,g,d',nisg(i),nisd(i) 

	p=p3(i);s=s1(i)      !根据中间压力 p 和1点的 s 计算2点的 h 
    call PSFLSH (p,s,z,t,D,Dl,Dv,x,y,q,e,h,cv,cp,w,ierr,herr)
	h2(i)=h;t2(i)=t

	t=tm(i)+0.1;p=p3(i)  !根据中间压力 p 和对应饱和温度 T+0.1 后计算3点的 h v3 s3（通过+0.1让状态点稍微位于过热区）
    call TPFLSH (t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
    h3(i)=h;s3(i)=s   ! 0.05882 = 1/17      v3(i)=0.05882/d;


	p=p5(i) ; s=s3(i)  !根据冷凝压力 p 和 3点的 s3 计算4点的 h4
    call PSFLSH (p,s,z,t,D,Dl,Dv,x,y,q,e,h,cv,cp,w,ierr,herr)
	h4(i)=h;t4(i)=t

    write (*,*),'t2 t4',t2(i),t4(i)

    t2(i)=t1(i)+(t2(i)-t1(i))/nisd(i)
	t4(i)=tm(i)+(t4(i)-tm(i))/nisg(i)

	write (*,*),'t2* t4*',t2(i),t4(i)

	t=t2(i);p=p3(i)         !根据2点 p T 计算
    call TPFLSH (t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
    h2(i)=h;s2(i)=s
	t=t4(i);p=p5(i)         !根据4点 p T 计算
    call TPFLSH (t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
    h4(i)=h;s4(i)=s

    ! *********以上 2 ,3, 4 ****************************

	t=tm(i)-0.1;p=p3(i)  !根据中间压力 p 和对应饱和温度 T-0.1 后计算7点的 h v7 （通过-0.1让状态点稍微位于过冷区）
	call TPFLSH (t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
	h7(i)=h;s7(i)=s

	t=t77(i);p=p3(i)  !根据中间压力 p 和对应饱和温度 T-0.1 后计算7点的 h v7 （通过-0.1让状态点稍微位于过冷区）
	call TPFLSH (t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
	h77(i)=h;s77(i)=s

	t=tl-0.1;p=p1(i)  !根据中间压力 p 和对应饱和温度 T-0.1 后计算9点
	call TPFLSH (t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
	h9(i)=h;s9(i)=s

	t=tl+0.1;p=p1(i)  !根据中间压力 p 和对应饱和温度 T-0.1 后计算9点
	call TPFLSH (t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
	h10(i)=h;s10(i)=s


    ! *********以上 7 ,9,10,7'(77)****************************

	x6(i)=(h5(i)-h7(i))/(h3(i)-h7(i))            !6点干度
	x8(i)=(h77(i)-h9(i))/(h10(i)-h9(i))        !8点干度
	write(*,*),'x干度6，8',x6(i),x8(i)

	s6(i)=x6(i)*s3(i)+(1-x6(i))*s7(i)
    s8(i)=x8(i)*s10(i)+(1-x8(i))*s9(i)

	if(x8(i)<0) then
    s8(i)=s9(i)
	end if

    h6(i)=h5(i);h8(i)=h77(i)
    !write(*,*),'ss',s6(i),s8(i)

    ! *********以上 6 ,8****************************

	t=ta;p=pa         
    call TPFLSH (t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
    ha=h;sa=s
    write (*,*), 'a环境25C',ha,sa

    ! *********以上环境 ha sa****************************

    write (*,*),'pl tl ph th',p1(i),tl,p5(i),th

	! *********基本****************************

    ex1(i)=h1(i)-ha-ta*(s1(i)-sa)
    ex1_sat(i)=h1_sat(i)-ha-ta*(s1_sat(i)-sa)
	ex2(i)=h2(i)-ha-ta*(s2(i)-sa)
	ex3(i)=h3(i)-ha-ta*(s3(i)-sa)
	ex4(i)=h4(i)-ha-ta*(s4(i)-sa)
	ex5(i)=h5(i)-ha-ta*(s5(i)-sa)
	ex6(i)=h6(i)-ha-ta*(s6(i)-sa)
	ex7(i)=h7(i)-ha-ta*(s7(i)-sa)
	ex77(i)=h77(i)-ha-ta*(s77(i)-sa)
	ex8(i)=h8(i)-ha-ta*(s8(i)-sa)

    write (*,*),'EX',ex1_sat(i),ex1(i),ex2(i),ex3(i),ex4(i),ex5(i),ex6(i),ex7(i),ex77(i),ex8(i)

    ! *********以上 Ex 计算****************************

    !p1=1000*p1; pm(i)=1000*p3(i); p5=1000*p5  
	! ;  h1=1.0/1000.0*h1; h2(i)=1.0/1000.0*h2(i)         !单位换算kPa-Pa，J/kg - kJ/kg
	!h3(i)=1.0/1000.0*h3(i); h4(i)=1.0/1000.0*h4(i); h5=1.0/1000.0*h5; h7(i)=1.0/1000.0*h7(i)  

    write (*,*), 'h焓',h1(i),h2(i),h3(i),h4(i),h5(i),h6(i),h7(i),h77(i),h8(i),h9(i)
    write (*,*), 's熵',s1(i),s2(i),s3(i),s4(i),s5(i),s6(i),s7(i),s77(i),s8(i),s9(i)
    write (*,*),'***********************************中间t和p', tm(i),p3(i)
    end do

   
    do i=1,n-1      !公式循环
	q0(i)=h1_sat(i)-h8(i)           !单位质量制冷量
	wd(i)=h2(i)-h1(i)           !低压压缩机比功

      !qvd(i)=qmd(i)*v1/0.65    !理论输气量 制冷原理书p126，p127的两个系数0.73和0.65          
	pd(i)=qmd(i)*wd(i)       !低压压缩机功率
    wg(i)=h4(i)-h3(i)        !高压压缩机比功
	qmg(i)=(h2(i)-h7(i))/(h3(i)-h5(i))   !高压压缩机进气质量流量 qmg*h6 + qmd*h2 = qmg*h3 + qmd*h7 (h5=h6)
	pg(i)=qmg(i)*wg(i)       !高压压缩机功率
	cop(i)=qmd(i)*q0(i)/(pd(i)+pg(i))  
	  !qvg(i)=qmg(i)*v3(i)/0.73

    mass_mol=capacity/(h1_sat(i)-h8(i))       !1p换算(此行以上的qmg(i)=单位1=1mol/s),以下为实际流量
    qmg(i)=qmg(i)*mass_mol
    qmd(i)=qmd(i)*mass_mol

 
	!xv(i)=qvg(i)/qvd(i)     !高低压进气体积比
	xm=qmg(i)/qmd(i)        !高低压进气质量流量比
	qk(i)=qmg(i)*(h4(i)-h5(i)) !冷凝器热负荷
	qeva(i)=qmd(i)*(h1_sat(i)-h8(i)) !蒸发器热负荷


    ! *********以下 Ex 损失计算****************************

    exloss_d(i)=qmd(i)*(wd(i)-(ex2(i)-ex1(i)))
	exloss_g(i)=qmg(i)*(wg(i)-(ex4(i)-ex3(i)))
    exloss_con(i)=qmg(i)*(ex4(i)-ex5(i))-qk(i)*(ex2_con-ex1_con)/(ha2_con-ha1_con)
    exloss_eva(i)=qmd(i)*(ex8(i)-ex1_sat(i))-qeva(i)*(ex2_eva-ex1_eva)/(ha2_eva-ha1_eva)
    exloss_thrg(i)=qmg(i)*ta*(s6(i)-s5(i))
    exloss_thrd(i)=qmd(i)*ta*(s8(i)-s77(i))
    exloss_sep(i)=qmd(i)*ex2(i)+qmg(i)*ex6(i)-qmg(i)*ex3(i)-qmd(i)*ex7(i)
	exloss_back(i)=qmd(i)*(ex7(i)+ex1_sat(i)-ex77(i)-ex1(i))

	ex_eff(i)=(qeva(i)*(ta/tl-1))/(qmg(i)*wg(i)+qmd(i)*wd(i))


	write (1,20),i,tm(i),p3(i)*1000,xm(i),cop(i)
	!write(*,*),'**************蒸发器用损失',ex8(i)-ex1(i),qeva(i)*(ta/tl-1)

    end do


10  format('序号   中间温度   中间压力   高低压质量流量比   COP   ')
20  format(I2,'      ',F6.2,'  ',F10.2,'        ',F6.4,'       ',F5.3)

	write (1,*), '以1p制冷量 2500 W为标准，计算各部件火用损失(W)' ; write (1,30)
30	format('序号 低压压缩机 高压压缩机   冷凝器    节流阀g   气液分离器   节流阀d    蒸发器    回热器    用效率')
40  format(I2,'    ',F7.3,'    ',F7.3,'    ',F7.3,'    ',F7.3,'    ',F7.3,'    ',F7.3,'    ',F7.3,'    ',F7.3,'    ',F7.3)   

    do i=1,n-1  
    write (1,40),i,exloss_d(i),exloss_g(i),exloss_con(i),exloss_thrg(i),exloss_sep(i),exloss_thrd(i),exloss_eva(i),exloss_back(i),ex_eff(i)
	!write (1,*),i,exloss_d(i),exloss_g(i),exloss_con(i),exloss_thrg(i),exloss_sep(i),exloss_thrd(i),exloss_eva(i)
	end do

    do i=1,n-1
	!write(*,*),'t1*  t1, t7*  t7',t1(i),tl,t77(i),tm(i)
	write(*,*),'x干度6，8',x6(i),x8(i)
	end do

	close(1)
	end program main

