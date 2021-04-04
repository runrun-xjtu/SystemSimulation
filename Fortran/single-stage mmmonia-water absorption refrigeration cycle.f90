!********************程序开始********************

	program main

    implicit double precision (a-h,o-z)
    implicit integer (i-k,m,n)
    parameter (ncmax=20)  
    dimension x(ncmax),xliq(ncmax),xvap(ncmax),f(ncmax)
    character hrf*3
    character*255 hf(ncmax),hfmix,hflnme,herr
    dimension xdew(ncmax),ybub(ncmax)
	
!********************声明部分********************

	integer i,j,n,iteration,judge(30)              ! 循环-判断 参量
    real temp(100)

	real e,Th1,Tc1,cph,cpc,qmh,qmc                     ! 换热器模块的输入
	real ntu,kA,Th2,Tc2,dtm                            ! 换热器模块的输出
	real qcool_con,qcool_abs,cpH2O,Tcool1,qmcool       ! heatx 2 冷却水模块(pcool2 与内部函数避免冲突 )
	real Qcapacity_origin                              ! heatx 3 & 蒸发器 模块
    real Tcool2                                        ! 吸收器 模块

    real mD,mW,mL,qcon,qgen,TD,TD_con,hIn,hW,hD_con,xW    ! 精馏塔模块的输出: 质量流量  各点焓值 塔顶温度(冷凝前后)
    real wF,TF,mF                                         ! 精馏塔模块的输入

	real wpump,pH,pL,pgen,pabs,h1,h2,h3,h4,h5,h6,h7,h8,h9,h911,h10,h11   ! 主程序所需变量
	real T2,T3,T4,T5,T6,T7,T8,T9,T911,T10,T11,T12   ! 主程序所需变量(T1为非自由设定变量)
    real cp2,cp3,cp4,cp5,cp6,cp8,cp9,cp10,m5,x5,x1,x2,COP   ! 主程序所需变量
   
	real :: T1=313.15,m1=10,a=10,pair=101.325      ! 非自由设定输入(a(相对挥发度)定值,被覆盖的假定输入)       m(质量流量) - g/s mol(摩尔流量) - mol/s                                                      
    real :: w1=0.3,TH=308.15,TL=253.15,Qcapacity=2500,wD=0.999,TW=403.15,Tcool0=298.15,e1=0.9,e2=0.9,e3=0.8,e4=0.9    ! 全循环 自由设定输入

!--------------------------------------------------
!********************主程序开始********************
!--------------------------------------------------
    open(1,file='result.txt')

  do iteration = 1,20

    !----- 冷凝T,p 蒸发T,p -----
	write(*,*),'----- ----- 冷凝T,p 蒸发T,p ----- -----'
	write(*,*)

    nc=1                     
    hf='ammonia.fld'      
    hfmix='hmx.bnc'          
    hrf='DEF'                 
    call SETUP (nc,hf,hfmix,hrf,ierr,herr)

	kph=2;t=TL
    call SATT (t,x,kph,p,rhol,rhov,xliq,xvap,ierr,herr);
    pL=p;
	kph=2;t=TH
    call SATT (t,x,kph,p,rhol,rhov,xliq,xvap,ierr,herr)
	pH=p;

	pgen=pH; pabs=pL
    write (*,*),'冷凝 T =',TH,'蒸发T=',TL,'(K)'
	write (*,*),'冷凝 p =',pH,'蒸发p=',pL,'(kPa)'

    !----- pump -----
	write(*,*)
	write(*,*),'----- ----- pump ----- -----'
	write(*,*)

    hflnme='am1wa.mix' ! 混合物:氨-水
    hfmix='hmx.bnc'
    hrf='DEF'
	call SETMIX (hflnme,hfmix,hrf,ncc,hf,x,ierr,herr)


    x1=(w1/17)/(w1/17+(1-w1)/18)
    x(1)= (w1/17)/(w1/17+(1-w1)/18) ; x(2)=1-x(1)
    t=T1; p=pabs
	call TPFLSH (t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
	h1=h ; h1=h1/17*x(1)+18*(1-x(1))
    t=T1; p=pgen
	call TPFLSH (t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
	h2=h  ; h2=h2/17*x(1)+18*(1-x(1))
	wpump=(h2-h1)*m1
	write (*,*),'NH3 摩尔分数(浓溶液)',x(1)
    write (*,*),'泵功(W)',wpump
	write (*,*),'泵前焓：',h1,'泵后焓：',h2,'(kJ/kg)'

	!----- rectifying tower & Heatx 1 -----
	write(*,*)
    write(*,*),'----- ----- rectifying tower & Heatx 1 ----- -----'
	write(*,*)

    wF=w1; TF=T1 ;mF=m1
	call rectification(wF,mF,pgen,TW,TF,wD,a,mD,mW,mL,qcon,qgen,TD,TD_con,hIn,hW,hD_con,xW)
    T5=TW; h5=hW ; m5=mW ;  x5=xW        ! 返回 T5,h5 用于溶液热交换器计算
	!write(*,*),T5,h5,m5,x5

	!----- Heatx 1 -----
	write(*,*)
	write(*,*),'----- ----- Heatx 1 ----- -----'
	write(*,*)

	x2=(w1/17)/(w1/17+(1-w1)/18)
	!write(*,*),3,3,x5,x2,e1,T5,T1,pgen,pgen,m5,m1
	call newhex(3,3,x5,x2,e1,T5,T1,pgen,pgen,m5,m1,ntu,kA,Th2,Tc2,dtm)  ! 调用 newhex 计算溶液热交换器
	wF=w1; TF=Tc2; T10=Th2; mF=m1 

    write(*,*)
	write(*,*),'----- ----- 重新调用 rectification ----- -----'
	write(*,*)
	call rectification(wF,mF,pgen,TW,TF,wD,a,mD,mW,mL,qcon,qgen,TD,TD_con,hIn,hW,hD_con,xW)

	!----- Heatx 2 -----
	write(*,*)
	write(*,*),'----- ----- Heatx 2 ----- -----'
	write(*,*)

    qcool_con=qcon

	nc=1                     
    hf='water.fld'      
    hfmix='hmx.bnc'          
    hrf='DEF'                 
    call SETUP (nc,hf,hfmix,hrf,ierr,herr)

    t=Tcool0 ; p=pair 
    call TPFLSH (t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
    cpH2O=cp ; cpH2O = cpH2O/18 ! J / mol K 转化为 J / g K

	Tcool1=Tcool0+e2*(TD_con-Tcool0)
	qmcool=(qcool_con/cpH2O)/(Tcool1-Tcool0)

	write(*,*),'冷凝器 冷却水温升(K):',(Tcool1-Tcool0),'冷却水流量(kg/s):',qmcool/1000
	write(*,*),'冷凝器冷却负荷(W):',qcool_con,'冷却水比热(J/(gK)):',cpH2O

	!----- Heatx 3 & evaporator -----
	write(*,*)
    write(*,*),'----- ----- Heatx 3 & evaporator ----- -----'
	write(*,*)

    nc=1                     
    hf='ammonia.fld'      
    hfmix='hmx.bnc'          
    hrf='DEF'                 
    call SETUP (nc,hf,hfmix,hrf,ierr,herr)

	call newhex(1,1,1.0,1.0,e3,TD_con,TL,pH,pL,mD,mD,ntu,kA,Th2,Tc2,dtm)  ! 调用 newhex 计算过冷器
	T9=Tc2 ; T6=Th2

    t=TL+0.1 ; p=pL
    call TPFLSH (t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
    h8=h; h8=h/18
	t=T6 ; p=pH
    call TPFLSH (t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
    h6=h; h6=h/18 

    Qcapacity_origin=mD*(h8-h6)
	COP=Qcapacity_origin/(qcon+wpump)

	write(*,*),'制冷量计算(初值):',Qcapacity_origin,'(W)'
	write(*,*),'COP计算(初值):',COP

	!----- absorber -----
	write(*,*)
    write(*,*),'----- ----- absorber ----- -----'
	write(*,*)

    t=T9 ; p=pL
    call TPFLSH (t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
    h9=h 

	hflnme='am1wa.mix' ! 混合物:氨-水  
    hfmix='hmx.bnc'
    hrf='DEF'
	call SETMIX (hflnme,hfmix,hrf,ncc,hf,x,ierr,herr)

    x(1)=xW ; x(2)=1-xW ; t=T10; p=pgen  !  计算10点焓
	call TPFLSH (t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
	h10=h 

	h911=(mW*h10+mD*h9)/mF

	x(1)=(w1/17)/(w1/17+(1-w1)/18) ; x(2)=1-x(1) ; p=pabs ; h=h911
    call PHFLSH (p,h,x,t,D,Dl,Dv,x,y,q,e,s,cv,cp,w,ierr,herr)
    T911=t


	call newhex(3,2,x1,1.0,e4,T911,Tcool1,pabs,pair,mF,qmcool,ntu,kA,Th2,Tc2,dtm)
	Tcool2=Tc2 ; T12=Th2
	
	write(*,*)
	write(*,*),'吸收器混合温度:',T911,'出口温度:',T12
	
	!----- iteration result-----
	write(*,*)
    write(*,*),'----- ----- iteration result ----- -----'
	write(*,*)

    if(abs(T12-T1)<0.01) then 
	   write(*,*),'循环已收敛 , |T12-T1|= ',abs(T12-T1)
       exit
	else 
       write(*,*),'T1=',T1,'T12=',T12,'|T12-T1|= ',abs(T12-T1)
	   pause   ! 中断检测
	   T1=T12
	end if

  end do

	

!--------------------------------------------------
!********************主程序结束********************
!--------------------------------------------------
    close(1)
	write(*,*)
	write(*,*),'!!! 程序结束 !!!'
	write(*,*)
    pause

!********************内部函数********************

    contains 

    !subroutine newhex(X_h,X_c,xhot,xcool,e0,Th1,Tc1,phot,pcool,qmh,qmc,ntu,kA,Th2,Tc2,dtm)
    !subroutine rectification(wF,mF,pgen,TW,TF,wD,a,mD,mW,mL,qcon,qgen,TD,TD_con,hIn,hW,hD_con,xW)

	subroutine rectification(wF,mF,pgen,TW,TF,wD,a,mD,mW,mL,qcon,qgen,TD,TD_con,hIn,hW,hD_con,xW)
	   real wF,mF,pgen,TW,TF,wD,a,mD,mW,mL,qcon,qgen,TD,TD_con,hIn,hW,hD_con,xW
       real xF,xD,wW,molF,molD,molW,molL  ! 精馏塔: F,D,W,L摩尔分数 质量分数, 流量(质量,摩尔)
	   real yq,yq_r,rmin,r         ! 精馏塔: q线,回流比(最小，操作)
       real hD,xy_zhu(101),jialiao,lilun           ! 精馏塔:各点焓值,逐板法参数,加料板数
	   real xNH3,tNH3(2),tpx,tin,pin

      hflnme='am1wa.mix' ! 混合物:氨-水
      hfmix='hmx.bnc'
      hrf='DEF'

    !-----xF,xD,xW,mD,mW (mF,wD,wF,pgen,TW)-----  
    
    xF=(wF/17)/(wF/17+(1-wF)/18)
	xD=(wD/17)/(wD/17+(1-wD)/18)

    call SETMIX (hflnme,hfmix,hrf,ncc,hf,x,ierr,herr)   ! 调用tpx函数(t,p->x)
	

    !-----嵌入tpx-----

     n=100;tin=TW;pin=pgen;tpx=-1

	   do i=1,n-1
          x(1)=1.0/n*i;x(2)=1-x(1);kph=1;p=pin
          call SATP (p,x,kph,t,rhol,rhov,xliq,xvap,ierr,herr);
          tNH3(1)=t
		  x(1)=1.0/n*(i+1);x(2)=1-x(1);kph=1;p=pin
          call SATP (p,x,kph,t,rhol,rhov,xliq,xvap,ierr,herr);
          tNH3(2)=t
		  if ( ((tin<tNH3(1)).AND.(tin>tNH3(2))) .OR. ((tin>tNH3(1)).AND.(tin<tNH3(2)))  ) then
		     tpx=1.0/n*(i+0.5)
			 write(*,10),tin,pin,tpx,i
		     exit
		  end if
       end do
10  format('tpx(t,p->x)调用完成，返回结果(t,p,x):',F7.3,' ',F8.3,' ',F7.3,' | ',I4,'次结束')

    !-----嵌入tpx-----
	xW=tpx
	
    wW=(xW*17)/(xW*17+(1-xW)*18)
	mD=mF*(wF-wW)/(wD-wW)
	mW=mF*(wF-wD)/(wW-wD)
    
    molF=mF/(17*xF+18*(1-xF))
	molD=mD/(17*xD+18*(1-xD))
	molW=mW/(17*xW+18*(1-xW))
	!temp(3)=mF*wF-mD*wD-mW*wW
	!temp(4)=molF*xF-molD*xD-molW*xW
	!write (*,*),temp(3),temp(4),molD,molW,molF
	write (*,*),'物料(g/s):D=',mD,'W=',mW,'xW=',xW

	!-----r,rmin-----

    yq_r=(a*xF)/(1+a*xF-xF)   ! 视为饱和液进料，q线与平衡线的交点(xF,yq_r)
    rmin=(xd-yq_r)/(yq_r-xF)
    r=1.5*rmin
	write (*,*),'回流比: rmin=',rmin,'r=',r

	!-----逐板计算法-----
    
	mL=r*mD;molL=r*molD

    write(*,*)
    write(*,*),'平衡线方程: y = ',a,' * x / ( 1 + ',a-1,' x )'
    write(*,*),'精馏段方程: y = ',r/(r+1),' * x + ',1.0/(r+1)*xD
    write(*,*),'提馏段方程: y = ',(molL+molF)/(molL+molD),' * x - ',(molW/(molL+molD))*xW
    write(*,*),'q线方程(饱和液进料): x = ',xF
	write(*,*),'精馏段 提馏段 交点坐标 ( ',xF,(r/(r+1))*xF+1.0/(r+1)*xD,')'
    write(*,*)
	!temp(1)=(r/(r+1))*xF+1.0/(r+1)*xD
	!temp(2)=((molL+molF)/(molL+molD))*xF-(molW/(molL+molD))*xW
	!temp(5)=molF-molD-molW
	!write (*,*),temp(1),temp(2),temp(5)
    
	judge(1)=0
    judge(2)=0
	xy_zhu(1)=xD

	do i=1,100                  ! 最多n+1个点,n次计算

	  if(mod(i,2)==1) then
        xy_zhu(i+1)=xy_zhu(i)/(a+xy_zhu(i)-a*xy_zhu(i))
		write(*,*),'第',i/2+1,'层板，计算得 x = ',xy_zhu(i+1)
        if ((xy_zhu(i+1)<xF).AND.(judge(2)==0)) then  ! judge(2)作用: 只运行一次
		  judge(1)=1 ; judge(2)=1
          jialiao=i/2+1
		  write(*,*),' x < xF 加料板:',jialiao
		end if
        if (xy_zhu(i+1)<xW) then
		  lilun=i/2+2
		  write(*,*),' x < xW 理论塔板数:(包括再沸器):',lilun
		  exit
		end if
	  end if

      if(mod(i,2)==0) then
	    if (judge(1)==0) then
	      xy_zhu(i+1)=(r/(r+1))*xy_zhu(i)+1.0/(r+1)*xD
		  !write(*,*),'y-jing',i,xy_zhu(i+1)
        end if
		if (judge(1)==1) then
          xy_zhu(i+1)=((molL+molF)/(molL+molD))*xy_zhu(i)-(molW/(molL+molD))*xW
		  !write(*,*),'y-ti',i,xy_zhu(i+1)
		end if
	  end if
	  
	end do

    !-----各点焓计算-----
    
	x(1)=xF; x(2)=1-xF; t=TF; p=pgen  ! 进料
	call TPFLSH (t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
	hIn=h

    x(1)=xW; x(2)=1-xW; t=TW; p=pgen  ! 塔釜
	call TPFLSH (t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
	hW=h

	x(1)=xD; x(2)=1-xD;kph=2;p=pgen  ! 塔顶 (冷凝前)
    call SATP (p,x,kph,t,rhol,rhov,xliq,xvap,ierr,herr);
    TD=t
    t=TD; p=pgen
	call TPFLSH (t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
	hD=h

      nc=1        ! 冷凝后视为纯NH3             
      hf='ammonia.fld'      
      hfmix='hmx.bnc'          
      hrf='DEF' 
	  call SETUP (nc,hf,hfmix,hrf,ierr,herr)  

	kph=1;p=pgen   ! 塔顶 (冷凝后)
	call SATP (p,x,kph,t,rhol,rhov,xliq,xvap,ierr,herr)
	TD_con=t
	t=TD_con-0.1; p=pgen
	call TPFLSH (t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
	hD_con=h
	
    write(*,*)  ! h - J/mol 变为 kJ/kg
	write(*,*),'进料温度:',TF,'塔顶温度(冷凝前):',TD,'(K)'
	write(*,*),'塔顶温度(冷凝后):',TD_con,'塔釜温度:',TW,'(K)'
	write(*,*),'进料焓(kJ/kg):',hIn/(17*xF+18*(1-xF)),'塔釜出料焓:',hW/(17*xW+18*(1-xW))
	write(*,*),'塔顶蒸汽焓(冷凝前)(kJ/kg):',hD/(17*xD+18*(1-xD)),'塔顶出料焓(冷凝后):',hD_con/(17*xD+18*(1-xD))
	  
    !-----热负荷计算-----

    qcon=molD*(hD-hD_con)*(1+r)   ! mol/s * J/mol = J/s = W
    qgen=qcon+molD*hD_con+molW*hW-molF*hIn
	write(*,*),'冷凝器热负荷(W):',qcon,'发生器热负荷(W):',qgen

    end subroutine rectification

!--------------------新--换热器模块-----------------------------
	subroutine newhex(X_h,X_c,xhot,xcool,e0,Th1,Tc1,phot,pcool,qmh,qmc,ntu,kA,Th2,Tc2,dtm)
	   real xhot,xcool,e0,Th1,Tc1,phot,pcool,qmh,qmc,ntu,kA,Th2,Tc2,dtm
	   real ch,cc,cmin,q1,c_ratio,cph1,cph2,cpc1,cpc2,cphot,cpcool,Maveh,Mavec,Mave
	   integer X_h,X_c,ii

       select case(X_h)   ! 1 - NH3 ; 2 - H2O ; 3 - NH3+H2O ; 4 - Air
	      case(1)

		    Maveh=17

            nc=1        ! 冷凝后视为纯NH3             
            hf='ammonia.fld'      
            hfmix='hmx.bnc'          
            hrf='DEF' 
	        call SETUP (nc,hf,hfmix,hrf,ierr,herr) 

		    t=Th1 ;p=phot 
            call TPFLSH (t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
            cph1=cp

			if(abs(cph1)>1000) then    ! 防 cp over undefined 机制
			   do ii=1,100
			      t=Th1-ii*0.5 ;p=phot
				  call TPFLSH (t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
                  cph1=cp
				  if(abs(cph1)<1000) then
				     write(*,*), '** cpover **'
				     exit
				  end if
			   end do
			end if

		  case(2)

		    Maveh=18

		  	nc=1                     
            hf='water.fld'      
            hfmix='hmx.bnc'          
            hrf='DEF'                 
            call SETUP (nc,hf,hfmix,hrf,ierr,herr)

		    t=Th1 ;p=phot 
            call TPFLSH (t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
            cph1=cp

			if(abs(cph1)>1000) then    ! 防 cp over undefined 机制
			   do ii=1,100
			      t=Th1-ii*0.5 ;p=phot
				  call TPFLSH (t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
                  cph1=cp
				  if(abs(cph1)<1000) then
				     write(*,*), '** cpover **'
				     exit
				  end if
			   end do
			end if

		  case(3)

		    Maveh=17*xhot+18*(1-xhot) 

  		    hflnme='am1wa.mix'
            hfmix='hmx.bnc'
            hrf='DEF'
            call SETMIX (hflnme,hfmix,hrf,ncc,hf,x,ierr,herr)

		    x(1)=xhot ; x(2)=1-xhot ; t=Th1 ;p=phot 
            call TPFLSH (t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
            cph1=cp

			if(abs(cph1)>1000) then    ! 防 cp over undefined 机制
			   do ii=1,100
			      t=Th1-ii*0.5 ;p=phot
				  call TPFLSH (t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
                  cph1=cp
				  if(abs(cph1)<1000) then
				     write(*,*), '** cpover **'
				     exit
				  end if
			   end do
			end if

          case(4)

		    Mave=28.96

		    nc=1                     
            hf='air.ppf'      
            hfmix='hmx.bnc'          
            hrf='DEF'                 
            call SETUP (nc,hf,hfmix,hrf,ierr,herr)

		    t=Th1 ;p=phot 
            call TPFLSH (t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
            cph1=cp

			if(abs(cph1)>1000) then    ! 防 cp over undefined 机制
			   do ii=1,100
			      t=Th1-ii*0.5 ;p=phot
				  call TPFLSH (t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
                  cph1=cp
				  if(abs(cph1)<1000) then
				     write(*,*), '** cpover **'
				     exit
				  end if
			   end do
			end if

	   end select

	   

!-----------------------------
	   select case(X_c)   ! 1 - NH3 ; 2 - H2O ; 3 - NH3+H2O ; 4 - Air
	      case(1)

		    Mavec=17

            nc=1        ! 冷凝后视为纯NH3             
            hf='ammonia.fld'      
            hfmix='hmx.bnc'          
            hrf='DEF' 
	        call SETUP (nc,hf,hfmix,hrf,ierr,herr) 

		    t=Tc1 ;p=pcool 
            call TPFLSH (t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
            cpc1=cp

		  case(2)

		    Mavec=18

		  	nc=1                     
            hf='water.fld'      
            hfmix='hmx.bnc'          
            hrf='DEF'                 
            call SETUP (nc,hf,hfmix,hrf,ierr,herr)

		    t=Tc1 ;p=pcool 
            call TPFLSH (t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
            cpc1=cp

		  case(3)

		    Mavec=17*xcool+18*(1-xcool)

  		    hflnme='am1wa.mix'
            hfmix='hmx.bnc'
            hrf='DEF'
            call SETMIX (hflnme,hfmix,hrf,ncc,hf,x,ierr,herr)

		    x(1)=xcool ; x(2)=1-xcool ; t=Tc1 ;p=pcool 
            call TPFLSH (t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
            cpc1=cp

          case(4)
		    
			Mavec=29

		    nc=1                     
            hf='air.ppf'      
            hfmix='hmx.bnc'          
            hrf='DEF'                 
            call SETUP (nc,hf,hfmix,hrf,ierr,herr)

            t=Tc1 ;p=pcool 
            call TPFLSH (t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
            cpc1=cp

       end select
!!!!!!-----------------------do 循环开始
       write(*,*),'********** 换热器迭代开始 **********'
       do i = 1,100
	   
	   if (i==1) then
	      cphot=cph1
		  cpcool=cpc1
	   else 
	      if ( (abs(cphot-0.5*(cph1+cph2))<=0.0001) .AND. (abs(cpcool-0.5*(cpc1+cpc2))<=0.0001) ) then
             write(*,*),'第',i-1,'次收敛,偏差:',cphot-0.5*(cph1+cph2),cpcool-0.5*(cpc1+cpc2)
			 write(*,*),'********** 换热器迭代完成 **********'
             exit
		  end if
          !pause   ! 中断 (检查使用)
		  write(*,*),'换热迭代第',i-1,'次'
		  write(*,*),'cp(h,h1,h2),偏差=',cphot,cph1,cph2,cphot-0.5*(cph1+cph2)
		  write(*,*),'cp(c,c1,c2),偏差=',cpcool,cpc1,cpc2,cpcool-0.5*(cpc1+cpc2)
	      cphot=0.5*(cph1+cph2)
		  cpcool=0.5*(cpc1+cpc2)

	   end if

       ch=qmh*cphot
	   cc=qmc*cpcool
!!!!!!!!!!!-----------------------------
	   if(ch>cc) then
	     cmin=cc
		 Tc2=Tc1+(Th1-Tc1)*e0
		 q1=qmc*cpcool*(Tc2-Tc1)
		 Mave=Mavec
		 Th2=Th1-q1/ch
	   else
	     cmin=ch
		 Th2=Th1-(Th1-Tc1)*e0
		 q1=qmh*cphot*(Th1-Th2)
		 Mave=Maveh
		 Tc2=Tc1+q1/cc
	   end if

!----------------------------- 第二次 select case
!----------------------------- 第二次 select case

       select case(X_h)   ! 1 - NH3 ; 2 - H2O ; 3 - NH3+H2O ; 4 - Air
	      case(1)

            nc=1        ! 冷凝后视为纯NH3             
            hf='ammonia.fld'      
            hfmix='hmx.bnc'          
            hrf='DEF' 
	        call SETUP (nc,hf,hfmix,hrf,ierr,herr) 

		    t=Th2 ;p=phot 
            call TPFLSH (t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
            cph2=cp

		  case(2)

		  	nc=1                     
            hf='water.fld'      
            hfmix='hmx.bnc'          
            hrf='DEF'                 
            call SETUP (nc,hf,hfmix,hrf,ierr,herr)

		    t=Th2 ;p=phot 
            call TPFLSH (t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
            cph2=cp

		  case(3)

  		    hflnme='am1wa.mix'
            hfmix='hmx.bnc'
            hrf='DEF'
            call SETMIX (hflnme,hfmix,hrf,ncc,hf,x,ierr,herr)

		    x(1)=xhot ; x(2)=1-xhot ; t=Th2 ;p=phot 
            call TPFLSH (t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
            cph2=cp

          case(4)

		    nc=1                     
            hf='air.ppf'      
            hfmix='hmx.bnc'          
            hrf='DEF'                 
            call SETUP (nc,hf,hfmix,hrf,ierr,herr)

		    t=Th2 ;p=phot 
            call TPFLSH (t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
            cph2=cp

		end select

!-----------------------------
	   select case(X_c)   ! 1 - NH3 ; 2 - H2O ; 3 - NH3+H2O ; 4 - Air
	      case(1)

            nc=1        ! 冷凝后视为纯NH3             
            hf='ammonia.fld'      
            hfmix='hmx.bnc'          
            hrf='DEF' 
	        call SETUP (nc,hf,hfmix,hrf,ierr,herr) 

		    t=Tc2 ;p=pcool 
            call TPFLSH (t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
            cpc2=cp

			if(abs(cpc2)>1000) then    ! 防 cp over undefined 机制
			   do ii=1,100
			      t=Tc2-ii*0.5 ;p=pcool
				  call TPFLSH (t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
                  cpc2=cp
				  if(abs(cpc2)<1000) then
				     write(*,*), '** cpover **'
				     exit
				  end if
			   end do
			end if


		  case(2)

		  	nc=1                     
            hf='water.fld'      
            hfmix='hmx.bnc'          
            hrf='DEF'                 
            call SETUP (nc,hf,hfmix,hrf,ierr,herr)

		    t=Tc2 ;p=pcool 
            call TPFLSH (t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
            cpc2=cp

			if(abs(cpc2)>1000) then    ! 防 cp over undefined 机制
			   do ii=1,100
			      t=Tc2-ii*0.5 ;p=pcool
				  call TPFLSH (t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
                  cpc2=cp
				  if(abs(cpc2)<1000) then
				     write(*,*), '** cpover **'
				     exit
				  end if
			   end do
			end if

		  case(3)

  		    hflnme='am1wa.mix'
            hfmix='hmx.bnc'
            hrf='DEF'
            call SETMIX (hflnme,hfmix,hrf,ncc,hf,x,ierr,herr)

		    x(1)=xcool ; x(2)=1-xcool ; t=Tc2 ;p=pcool 
            call TPFLSH (t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
            cpc2=cp

			if(abs(cpc2)>1000) then    ! 防 cp over undefined 机制
			   do ii=1,100
			      t=Tc2-ii*0.5 ;p=pcool
				  call TPFLSH (t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
                  cpc2=cp
				  if(abs(cpc2)<1000) then
				     write(*,*), '** cpover **'
				     exit
				  end if
			   end do
			end if

          case(4)

		    nc=1                     
            hf='air.ppf'      
            hfmix='hmx.bnc'          
            hrf='DEF'                 
            call SETUP (nc,hf,hfmix,hrf,ierr,herr)

            t=Tc2 ;p=pcool 
            call TPFLSH (t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
            cpc2=cp

			if(abs(cpc2)>1000) then    ! 防 cp over undefined 机制
			   do ii=1,100
			      t=Tc2-ii*0.5 ;p=pcool
				  call TPFLSH (t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
                  cpc2=cp
				  if(abs(cpc2)<1000) then
				     write(*,*), '** cpover **'
				     exit
				  end if
			   end do
			end if

       end select

	   end do

!-----------------------------do 循环结束
 
	   dtm=((Th1-Tc2)-(Th2-Tc1))/alog((Th1-Tc2)/(Th2-Tc1))
	   !kA=q/dtm
	   !ntu=kA/cmin
       
	   if(ch>cc) then
	   c_ratio = cc/ch
	   else 
       c_ratio = ch/cc
	   end if
	   ntu=(-1)*alog((1-e0)/(1-e0*c_ratio))/(1-c_ratio)
	   kA=cmin*ntu

	   write (*,*),'热流体进口(K):',Th1,'热流体出口(K):',Th2
	   write (*,*),'冷流体进口(K):',Tc1,'冷流体出口(K):',Tc2
	   write (*,*),'对数平均温差(K):',dtm,'传热量(W):',q1/Mave
	   write (*,*),'NTU:',ntu,'KA (J/(Ks)):',kA/Mave

    end subroutine newhex
          
!********************程序结束********************

	end program main
