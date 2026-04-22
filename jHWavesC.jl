"""
	jHWavesC

This set of functions computes coastal-trapped wave properties in a stratified ocean.
This collection seeks a complex frequency.

To use this collection:

Make an input tuple:

	arr = jHWavesCSetup()

Modify an input tuple:

	narr = jHWavesCFinch(arr)

To run:

	jHWavesCM(arr)


Translated from Matlab bigc___ code  (6/8/2018 version)

K.H. Brink  4/20/2026
"""



"""
	Function jHWavesCM(arr)

	Main function for computing stratified coastal-trapped wave modes
		with complex frequency.

	call as
		jHWavesCM(arr)

	When it is done, it will ask if you want to save the results

K.H. Brink 4/20/2026

"""



function jHWavesCM(arr)

#   Julia version file to compute coastal trapped wave modes 
#       with stratification, bottom topography, and mean flow.
#       Allows for complex frequency.
#       Modes are not usually orthogonal.       
#       Calls various other functions beginning with "jHWavesC" or "jHWaves"
#
#   Call as 
#       jHWavesCM(arr)
#   where
#       arr is an input tuple containing all the pertinent data.
#
#   See function jHWavesCSetup to help set this up, and jHWavesCFinch to change the input
#       tuple conveniently
#
#   K. Brink  4/20/2026


# arr is a tuple consisting of:
#
#	nn	    number of x grid points
#	mm	    number of vertical grid points
#       wz          where wz is the first guess frequency (rad/sec), entered as an array
#       del         del = 0 for rigid lid; del = 1 for free surface
#       icbc        icbc = 0 says no flow through coastal wall at x = 0
#                   icbc !=0 says du/dx = 0 at shallow water boundary (requires flat bottom)
#       iobc        iobc = 0 says no flow through offshore boundary
#                   iobc = 1 says open offshore boundary
#       f           Coriolis parameter (rad/sec)
#       xmax        offshore size of grid  (km)
#       eps         nominal fractional accuracy for solutions (0.001 works)
#       npts        number of points on dispersion curve
#       rlz         first alongshore wavenumber (rad/cm)
#       drl         increment of alongshore wavenumber (rad/cm)
#       ndep        number of offshore distance, depth pairs to read  (in km, m)
#                       Must be >= 1
#       xdep   (ndep times)    offshore distances (km) for topography. First value must = 0
#       depr   (ndep times)    depth (m) corresponding to x
#       nr          Number of distance, bottom friction pairs (>=0)
#                       nr = 0 means no bottom friction
#       xr  (nr times)  offshore distance for r in km  
#       rr  (nr times)  bottom resistance coefficient values  in cm/sec
#       nnsq        number of Nsquared values to read  (>=1)
#       zr          depth increment for Nsquared to read  (m). First depth is
#                       at the surface.
#       alph        exponential tail scale (km) for extrapolating Nsquared beyond
#                       the max depth read
#       nsqr    (nnsq times)    Nsquared values (rad^2/s^2)
#       vzero       amplitude of mean alongshore flow (cm/sec)
#                      If vzero = 0, mean flow is zero everywhere
#       xzero       offshore distance of maximum in alongshore mean flow (km)
#       zzero       depth of mean flow maximimum  (m) (> 0 for subsurface max)
#       zscaled     downward e-folding scale for mean velocity (m)
#       zscaleu     upward e-folding scale for mean velocity (m)
#       xscaleoff   offshore e-folding scale for mean velocity (km)
#       xscaleon    onshore e-folding scale for alongshore velocity (km)
#       kk          determines where Nsquared is undisturbed                  
#        		kk = 1          undisturbed offshore
#                       kk = 0          undisturbed at coast
#       ipause      pause to see graphics (1) or execute without pauses (0)


#       All internal works are in cgs units, although inputs are in convenient units

#       Numerical grid is arranged so there is a row of grid points outside the physical
#           domain on all boundaries. Thus, there are mm-2 vertical grid points within,
#           or on the border of, the physical domain.
#       The grid points are arranged so that it is preferable to have nn > mm



	global nn, mm, dt, dx, f, rl, del, icbc, iobc, BB, h


	#   Extract array size
	nn = arr[1]
	mm = arr[2]
	println("nn, mm = ", nn, ",  ", mm)

	#   First frequency guess
	wgr = arr[3]
	wg = wgr[1] + im*wgr[2]

	#   Rigid lid/free surface, coastal wall or open, etc.
	del = arr[4]
	icbc = arr[5]
	iobc = arr[6] 
	println(' ')

	if del > 0.5
    		println("Free surface")
	else
    		println("Rigid lid")
	end
	if icbc > 0.5
    		println("Open BC at x = 0")
	else
    		println("Closed BC at x = 0")
	end
	if iobc > 0.5
    		println("Open BC at x = xmax")
	else
	        println("Closed BC at x = xmax")
	end

	#   f
	f = arr[8]
	println(' ')
	@printf("f = %.3e rad/sec \n", f)

	#   Domain size offshore
	xmax = arr[9]

	#   Desired fractional accuracy
	eps = arr[10]

	#   number of dispersion curve points, first wavenumber and wavenumber increment
	npts = arr[11]
	rlz = arr[12]
	drl = arr[13]
	rll = rlz*ones(npts) + drl*(0:(npts-1))
	wnn = 0*ones(npts,2)

	#   Extract depth information (#, x, depth)
	ndep = arr[14]
	xdep = arr[15]
	depr = arr[16]

	#   Extract bottom friction information (number, x, r)
	nr  = arr[17]
    	xr = arr[18]
    	rr = arr[19]

	#   Extract base-state Nsquared information
	nnsq = arr[20]
	zr = arr[21]
	alph = arr[22]
	nsqr = arr[23]

	#   Pull out information on mean alongshore flow
	vzero = arr[24]
    	xzero = arr[25]
    	zzero   = arr[26]
    	zscaled = arr[27]
    	zscaleu = arr[28]
    	xscaleoff = arr[29]
    	xscaleon = arr[30]
    	kk = arr[31]

	#    Get the command for pausing during execution or not
	ipause = arr[32]

#   Now get started: calculate the mean fields

	xmax = xmax*1.0e05
	dt = 1/(mm-3)
	dx = xmax/(nn-3)


	icoastflag = jHWavesDep(xdep,depr)
	jHWavesConChk()
	jHWavesRr(nr,xr,rr)
	jHWavesNsq(zr,alph,nsqr) 

	vtup  = jHWavesVzCal(vzero,xzero,zzero,zscaleu,zscaled,xscaleoff,xscaleon)
	rtup = jHWavesRhoCal(vzero,vtup, kk,ipause)

	
	igravflag = rtup[2]


#   Check to see if there are flags for gravitational instability or inconsistency
	if igravflag !=0
		println(' ')
    		println("Exit because density is gravitationally unstable")
		return NaN
	end

	if icoastflag > 0.05 && icoastflag < 1.5
    		println(' ')	
    		println("Exit because of attempt to use open x = 0 boundary condition when the bottom is not flat there")	
		return NaN
	end

	if icoastflag > 1.5
		println(' ')
    		println("Exit because of attempt to use open x = xmax boundary condition when the bottom is not flat there")
		return NaN
	end

	sleep(1)

#   Set parameters for search
	maxit = 80
	wnn = NaN*ones(npts,1)*(1 + im)
	

	xgr = 0*ones(nn-2,mm-2)
    	zgr = 0*ones(nn-2,mm-2)
    	for n = 1:nn-2
        	xtemp = dx*(n-1)
        	xgr[n,:] = xtemp*ones(1,mm-2)/1e5
        	zgr[n,:] = h[n]*(-1*ones(mm-2) + dt*(0:mm-3))/100
    	end
	p = NaN*ones(nn-2,mm-2)*(1+im)
	BB = NaN*ones(nn*mm)*(1 + im)
	exitflag = NaN


#   Now actually do the calculations

	ncal = 1

	while ncal <= npts
		GLMakie.closeall()
    		rl = rll[ncal]
    		println(' ')
    		@printf("Wavenumber = %.3e rad/cm", rl)
    		println(' ')
		println(' ')
 		println("Iterations:")  

		ggg = [real(wg), imag(wg)]
		steppr = 0.025*abs.(real(wg))						#	search simplex
		steppi = 0.025*abs.(imag(wg))	
		ds = [steppr,  steppi]

		#steppi = 0.025*abs(wg)
		opt = Opt(:LN_NELDERMEAD,2)
		NLopt.min_objective!(opt,jHWavesCDrvr)
		NLopt.xtol_rel!(opt,eps)
		NLopt.maxeval!(opt,maxit)
		NLopt.initial_step!(opt, ds)


		#ggg = [wg]
		(minf, minx, ret) = NLopt.optimize!(opt,ggg)
		
		neval = NLopt.numevals(opt)

		if ret != :XTOL_REACHED								# ret is a "symbol" variable	
			exitflag = 0
			println(' ')
			println("Used up ", maxit, " iterations")
			www = NaN*(1 + im)
		else
			exitflag = 1
			www = minx[1] + minx[2]*im
		end


    
   # output result
   
   	 	if exitflag == 1
        		println(' ')
        		println("Converged!")
        		println(' ')
        		wnn[ncal,1] = www
     
        		if ncal ==1 && npts != 1
            			c = www/rl
            			wg = c*rll[2]
        		elseif ncal >=2
            			cg = wnn[ncal,1] - wnn[ncal-1,1]
            			wg = www + cg'
            			if ncal > 2.5
                			dc = wnn[ncal,1] - 2*wnn[ncal-1,1] + wnn[ncal-2,1]
                			wg = wg + dc
            			end
        		end
        		

			fig3 = Figure()
			ax1 = Axis(fig3[1,1])
			lines!(ax1,rll,real(wnn[:,1]), color = :black)
			lines!(ax1,rll,imag(wnn[:,1]), color = :red)
			ax1.ylabel = "ω [1/s]"
			ax1.xlabel = "Wavenumber [1/cm]"
			ax1.title = "Dispersion Curve (Imaginary in red)"
			xlims!(ax1,0.,maximum(rll))
			hidedecorations!(ax1; label = false, ticklabels = false, ticks = false )
			display(fig3)


#		Calculate and plot pressure
        
       		
        		wj = jHWavesCDrvr(minx,0.)
			csign = conj(BB[(2*mm-1)])/abs2(BB[(2*mm -1)])
        		
    			BB = BB*csign

    			for n = 1:nn-2
        			mlow = n*mm + 2
        			mhigh = (n+1)*mm - 1
        			p[n,:] = BB[mlow:mhigh]
    			end
			p = jHWavesCNorm(p)


			if ipause > 0.5				

				xff = xgr[:,1]
				hh = -zgr[:,1]
				ybb = -maximum(h)*ones(nn-2)
    
				fig4 = Figure()
				ax1 = Axis(fig4[1,1])

				tol = 0.001
				rmm = maximum(maximum(abs.(real(p))))
				amm = maximum(maximum(abs.(imag(p))))
				rata = amm/rmm
				pfac = 1
				if rata <= tol
					pfac = NaN
				end

				civ = jHWavesCCon(p,10)

				contour!(ax1, xgr,zgr,real(p), levels =civ, color = :black)	
				contour!(ax1,xgr,zgr,real(p),levels = vec([0 0]), linewidth = 3, color = :black)
				contour!(ax1, xgr,zgr,pfac*imag(p), levels = civ, color = :red)	
				contour!(ax1,xgr,zgr,pfac*imag(p),levels = vec([0 0]), linewidth = 3, color = :red)


				band!(ax1,xff,ybb,-hh,color = :blue)
				xlims!(ax1,0,maximum(xgr[:,1]))
				ylims!(ax1,minimum(zgr[:,1]),0)
				
				ax1.xlabel = "x [km]"
				ax1.ylabel = "z [m]"
				ax1.title = "Pressure (heavy contour is 0)"
				hidedecorations!(ax1; label = false, ticklabels = false, ticks = false )
				
				display(GLMakie.Screen(),fig4)

				if ncal < npts
					println(" ")
					println("Sleeping for 10 seconds")
					sleep(10)
				end	
			end
			ncal = ncal + 1
			
		else
			ncal = npts + 1
			
		end
	end
									# end of ncal loop

	rtup = jHWavesRhoCal(vzero,vtup, kk,0)



	
	if ipause  < 0.5 && exitflag == 1
		p = jHWavesCNorm(p)

		xff = xgr[:,1]
		hh = -zgr[:,1]
		ybb = -maximum(h)*ones(nn-2)
    
		fig4 = Figure()
		ax1 = Axis(fig4[1,1])

		tol = 0.001
		rmm = maximum(maximum(abs.(real(p))))
		amm = maximum(maximum(abs.(imag(p))))
		rata = amm/rmm
		pfac = 1
		if rata <= tol
			pfac = NaN
		end


		contour!(ax1, xgr,zgr,real(p), levels = 10, color = :black)	
		contour!(ax1,xgr,zgr,real(p),levels = vec([0 0]), linewidth = 3, color = :black)
		contour!(ax1, xgr,zgr,pfac*imag(p), levels = 10, color = :red)	
		contour!(ax1,xgr,zgr,pfac*imag(p),levels = vec([0 0]), linewidth = 3, color = :red)

		band!(ax1,xff,ybb,-hh,color = :blue)
		xlims!(ax1,0,maximum(xgr[:,1]))
		ylims!(ax1,minimum(zgr[:,1]),0)
				
		ax1.xlabel = "x [km]"
		ax1.ylabel = "z [m]"
		ax1.title = "Pressure (heavy contour is 0)"
		hidedecorations!(ax1; label = false, ticklabels = false, ticks = false )
				
		display(GLMakie.Screen(),fig4)
	end			

	println(' ')
	if exitflag  == 1
		jHWavesCEngDiag(wnn[ncal-1,1],BB)
	end


#	Plot final dispersion curve (real and complex)
	
	if npts > 1.5 
		fig2 = Figure()
		ax1 = Axis(fig2[1,1])
		#ax2 = Axis(fig2[2,1])
		lines!(ax1,rll,vec(real.(wnn)), color = :black)
		lines!(ax1,rll,vec(imag.(wnn)), color = :red)
		ax1.ylabel = "ω [1/s]"
		ax1.xlabel = "Wavenumber [1/cm]"
		ax1.title = "Dispersion Curve"
		xlims!(ax1,0.,maximum(rll))
		hidedecorations!(ax1; label = false, ticklabels = false, ticks = false )

		display(fig2)
	
	end

    
#       Option for saving results

	isave = inputi("Do you want to save these results? 0 (no) or 1 (yes) ")
    	if isave > 0.5
		println("Enter the name of the desired file (.jld2 will be appended to this) ")	
		str = readline()
		fname = str*".jld2"
		ilong = NaN
		otup = [f, ilong, icbc, iobc, del, rll, wnn, xgr, zgr, p, xr, rr, vzero]
		save_object(fname,otup)
		
		println(' ')
		println("Saved as a tuple in file ", fname)
		verbiage = """
			Contents of the tuple:

				1) f	(Coriolis parameter: 1/sec)
				2) 	(empty)
				3) icbc	 (= 0 for coastal wall, = 1 for open x = 0 BC)
				4) iobc	 (= 0 for wall at offshore boundary, = 1 for open)
				5) del	 (= 0 for rigid lid, = 1 for free surface)
				6) rll (wavenumber)
				7) wnn (complex frequency)
				8) xgr (x grid in km)
				9) zgr (z grid in m)
				10) p  (pressure)
				11) xr (x for bottom frictional coefficient: km)
				12) rr (bottom frictional coefficient: cm/sec)
				13) vzero  (extreme value of mean velocity: cm/sec)

			You can load this tuple with the 'load_object' command.
			"""
			
		println(verbiage) 
         
			
    	end

end				#  end of function




"""	
	jHWavesDep

	Calculate depth and its derivative

K.H. Brink 11/21/2025
"""


function jHWavesDep(xr,depr)

#       Compute the depth profile, given the inputs
#       Only compute depth for points within the domain or on the edges
#       iflag != 0 signals trouble with using open coastal or offshore boundary
#           i.e. hx != 0    
#
#       K.H. Brink 11/21/2025  based on  Matlab version of 5/2/03

	global nn, dx, h, hx, hxx, icbc, iobc


	x = dx*(0:(nn-3))
	xmax = x[nn-2]

	xr = xr*1.0e5
	depr = depr*100.

	if xmax > xr[end]
		xr[end] = xmax
	end

	itp = interpolate((vec(xr),),vec(depr),Gridded(Linear()))
	h = itp(x)

    	dxs = dx*dx
    	dx2 = 2*dx
	hx = 0*ones(nn-2)
	hxx = 0*ones(nn-2)

    	hx[2:(nn-3)] = (h[3:(nn-2)]-h[1:(nn-4)])/dx2
    	hxx[2:(nn-3)] = (h[3:(nn-2)] - 2*h[2:(nn-3)] + h[1:(nn-4)])/dxs

#       Assume the bottom is flat for x < 0, x > xmax

    	hx[1] = (h[2]-h[1])/dx2
    	hxx[1] = (-h[1] + h[2])/dxs
    
    	hx[nn-2] = (h[nn-2]-h[nn-3])/dx2
   	hxx[nn-2] = (-h[nn-3] +h[nn-2])/dxs
    
	iflag = 0
	if icbc != 0
    		if hx[1] != 0
        		iflag = 1
    		end
	end

	if iobc != 0
    		if hx[nn-2] != 0
        		iflag = 2
    		end
	end
	return iflag
end									# end of function





"""
	jHWavesConChk


	Check for consistency of bottom slope

K.H. Brink  11/21/2025
"""

function jHWavesConChk()

#       Check to find max of (hx/h)*(dx/dt)
# 
#    K.H. Brink 11/21/2025 based on Matlab version of  3/8/2004

	global hx, h, dx, dt

	arrt = (hx./(h + 0.5*hx*dx))*dx/dt
	arrt = abs.(arrt)

	rmax, ii = findmax(vec(arrt))

	println(' ')
	xm = dx*(ii-1)/1e5
	@printf(" Max consistency ratio = %.3f at x = %.2f km \n", rmax, xm)
	println("         This should be kept less than one, and definitely less than 10")
	println(' ')
end										# end of function





"""
	jHWavesRr(nr,xr,rr)

	Interpolate bottom friction coefficient onto the grid


K.H. Brink		11/21/2025
"""

function jHWavesRr(nr,xr,rr)

#   Compute the array of bottom friction parameters, given the inputs
#
#	K.H. Brink 11/21/2025, based on Matlab version of 5/5/2003

	global nn, dx, r, rx

	if nr == 0
      		r = 0*ones(nn-2)
      		rx = r
  	elseif nr ==1
      		r = rr[1]*ones(nn-2)
      		rx = 0*ones(nn-2)
  	else

    		x = dx*(0:(nn-3))
    		xmax = x[nn-2]
    		xr = xr*1.0e5
		r = 0*ones(nn-2)
      		rx = 0*ones(nn-2)

    		xrmax = xr[end]

    		ii = findall( x-> x > xrmax, x)
    		r[ii] = rr[end]*ones(length(ii))

    		ii = findall(x -> x <= xrmax, x)

		itp = interpolate((vec(xr),),vec(rr),Gridded(Linear()))
		r[ii] = itp(x[ii])

    		tdx = 2*dx

    		rx[2:(nn-3)] = (r[3:nn-2]- r[1:(nn-4)])/tdx

    		rx[nn-2] = 0
    		rx[1] = (r[2]-r[1])/dx

	end
end									# end of function







"""
	jHWavesNsq(zr,alph, nsqr)

	Compute undisturbed N^2 profile


K.H. Brink 11/21/2025
"""




function jHWavesNsq(zr,alph,nsqr)

#  Compute the undisturbed (no mean flow) N^2(z) profile
#       given the inputs
#
#   Gives back nsq(mm-2), the array of interior N^2 points, starting at the bottom
#       (consistent with grid point numbering)
#
#   K.H. Brink 11/21/2025 based on Matlab 5/5/2019 version

	global  mm, dt, h, nsq

	zr = zr*100
	alph = alph*1e5
	n = length(nsqr)

	nsq = ones(mm-2)
    	zmax = maximum(h)
    	nnmax = Int(ceil(zmax/zr) +2)
   	nsqtemp = 0*ones(nnmax)
    	zf = -zr*(0:(nnmax-1))
    	zrr = -zr*(0:(n-1))

    	if zmax > maximum(-zrr)
        	nsqtemp[1:n] = nsqr
        	ii = findall(zf -> zf < minimum(zrr),zf)
        	zfill = zf[ii]
        	nsqtemp[(n+1):nnmax] = nsqr[n]*exp.((zfill .- minimum(zrr))/alph)
        	zzz = zf
    	else
        	zzz = zrr
        	nsqtemp = nsqr
    	end
    	zint = -dt*zmax*(0:(mm-3))

	zzzt = reverse(zzz)
	nsqtempt = reverse(nsqtemp)
    
	itp = interpolate((vec(zzzt),),vec(nsqtempt), Gridded(Linear()))
	nsq = itp(zint)
   
    	nsq = nsq[(mm-2):-1:1]

     	ij = findall(nsq -> 0*nsq !=0, nsq)
    	if isempty(ij) == 0
        	zzzz = min(nsq)
        	nsq[ij] = zzzz
    	end
end							# end of function





""" 
	jHWavesVzCal(vzero,xzero,zzero,zscaleu,zscaled,xscaleoff,xscaleon)

	Compute the mean velocity field

K.H. Brink  11/21/2025
"""


function  jHWavesVzCal(vzero,xzero,zzero,zscaleu,zscaled,xscaleoff,xscaleon)

#   compute mean velocity field, given inputs
#   vmean,vmx,vmxx are in sigma coordinates
#   vtemp, vztemp in z coordinates
#       vtemp, vztemp to have size(nn-2,mm-2)    
# 
#       K.H. Brink 11/21/2015 based on Matlab version of 9/3/04

	global nn, mm, dt, dx, h, vm, vmx, vmxx

	xzero = xzero*1.0e5
	xscaleoff = xscaleoff*1.0e5
	xscaleon = xscaleon*1.0e5
	zzero = -zzero*100
	zscaleu = zscaleu*100
	zscaled = zscaled*100

	vm = 0*ones(nn-2,mm-2)
	vmx = 0*ones(nn-2,mm-2)
	vmxx = 0*ones(nn-2,mm-2)
	vtemp = 0*ones(nn-2,mm-2)
	vztemp = 0*ones(nn-2,mm-2)

	if abs(vzero) > 0.0001
#       Compute vmean on rectanguar grid first (vtemp, vztemp), then
#           interpolate into sigma (vm)

    		xfactoff = xscaleoff^2
    		xfacton = xscaleon^2
    		zfactu = zscaleu^2
    		zfactd = zscaled^2
    
    		hmax = maximum(h)
    		vtemp = NaN*ones(nn-2,mm-2)
    		x = dx*(0:(nn-3))
    		z = hmax*(-1*ones(mm-2) + dt*(0:(mm-3)))
    		dz = hmax*dt
    
    		ii = findall(x ->  x < xzero, x)
    		iic = findall(x -> x >= xzero, x)
    		jj = findall(z -> z < zzero, z)
    		jjc = findall(z -> z >= zzero, z)

    
   		if isempty(ii) != true
        		for n = 1:maximum(ii)
            			if isempty(jj) != true
					nzz = length(jj)
					zzn = zzero*ones(nzz)
                			vtemp[n,jj] = exp.(-((x[n]-xzero).*(x[n]-xzero)/xfacton))*exp.(-(z[jj]-zzn).*(z[jj]-zzn)/zfactd)
            			end
            
            			if isempty(jjc) != true
					nzz = length(jjc)
					zzn = zzero*ones(nzz)
                			vtemp[n,jjc] = exp.(-((x[n]-xzero).*(x[n]-xzero)/xfacton))*exp.(-(z[jjc]-zzn).*(z[jjc]-zzn)/zfactu)
            			end
        		end    
    		end
        
    		if isempty(iic) != true
        		for n = minimum(iic):nn-2
            			if isempty(jj) != true
					nzz = length(jj)
					zzn = zzero*ones(nzz)
                 			vtemp[n,jj] = exp.(-((x[n]-xzero).*(x[n]-xzero)/xfactoff))*exp.(-(z[jj]-zzn).*(z[jj]-zzn)/zfactd)
            			end
            
            			if isempty(jjc) != true
					nzz = length(jjc)
					zzn = zzero*ones(nzz)
                	 		vtemp[n,jjc] = exp.(-((x[n]-xzero).*(x[n]-xzero)/xfactoff))*exp.(-(z[jjc]-zzn).*(z[jjc]-zzn)/zfactu)
            			end          
        		end
        
    
    		end
		vtemp = vzero*vtemp
    
    		dz2 = 2*dz    	
		dx2 = 2*dx
    		dxs = dx*dx
    
    		vztemp[:,2:(mm-3)] = (vtemp[:,3:mm-2]-vtemp[:,1:(mm-4)])/dz2
    		vztemp[:,1] = (vtemp[:,2]-vtemp[:,1])/dz
    		vztemp[:,(mm-2)] = (vtemp[:,(mm-2)]-vtemp[:,(mm-3)])/dz

  
    		vmc = vtemp
    		vmxc = NaN*ones(nn-2,mm-2)
    		vmxxc = NaN*ones(nn-2,mm-2)
    
    		for n = 2:nn-3
        		vmxc[n,:] = (vmc[(n+1),:] - vmc[(n-1),:])/dx2
        		vmxxc[n,:] = (vmc[(n+1),:] - 2*vmc[n,:] + vmc[(n-1),:])/dxs
   		end
    
    		vmxc[1,:] = (vmc[2,:]-vmc[1,:])/dx
    		vmxxc[1,:] = 0*vmc[1,1:mm-2]
    
    		vmxc[nn-2,:] = (vmc[nn-2,:] - vmc[nn-3,:])/dx
    		vmxxc[nn-2,:] = 0*vmc[nn-2,1:mm-2]
   	
       

    		vm = jHWavesZtoSig(vmc)
    		vmx = jHWavesZtoSig(vmxc)
    		vmxx = jHWavesZtoSig(vmxxc)  

	end 

	vtup = (vtemp, vztemp)
	return vtup

end						# end of function







"""
	jHWavesRhoCal(vzero, vtup, kk, ipause)

	Calculate the density field etc., accounting for the mean flow



K.H. Brink 11/24/2025
"""

function jHWavesRhoCal(vzero,vtup,kk,ipause)

#   Compute the overall density field, given the basic N^2 and mean v
#   Then compute N^2, M^2 etc.
#	call as
#		rtup = jHWavesRhoCal(vzero, vtup, kk, ipause)
#
#   Does calculations in z coordinates and then converts to sigma
#
#   Arrays such as n2 have dimension (nn-2,mm-2)       
#  
#	K.H. Brink, 12/10/2025 based on Matlab of 5/21/2018

	global nn, mm, dt, dx, f, h, hx, n2, n2z, m2, m2z, m2x, nsq, vmx

	rho = 0*ones((nn-2),(mm-2))

	vtemp = vtup[1]
	vztemp = vtup[2]

	conts = -1.03/980
	cc = -f*1.03/980
	hmax = maximum(h)
	xmax = dx*(nn-3)/1e5		
	rhsave = 0.
	dz = hmax*dt

#   Compute background density  (as if vm = 0)
	rhob = 0*ones(mm-2)
	rhob[1:mm-2] = dz*conts*(cumsum(nsq[1:(mm-2)]) - 0.5*(nsq[1]*ones(mm-2) + nsq[1:(mm-2)]))
	rhob = rhob - minimum(rhob)*ones(mm-2)
	rhop = 0*ones(nn-2,mm-2)

	for n = 1:nn-2
    		rho[n,:] = rhob
	end

	if abs(vzero) > 0.001
    		for m = 1:(mm-2)
        		rhop[1,m] = 0.       
        		rhop[2:(nn-2),m] = cc*dx*(cumsum(vztemp[2:(nn-2),m]) + 0.5*( vztemp[1,m]*ones(nn-3) - vztemp[2:nn-2,m]))
        		if m == 1
            			(yyy,ii) = findmax(vec(abs.(rhop[:,m])))
            			rhsave = rhop[ii,m]
        		end
        		rhop[:,m] = rhop[:,m] - kk*rhop[nn-2,m]*ones(nn-2)
        		rho[:,m] = rho[:,m] + rhop[:,m]
    		end


    #   This assures that the maximum density change is not at the very bottomi
    		rho = rho - rhsave*ones((nn-2),(mm-2))
 
	
	end



	ccc = -980/1.03
	dz2 = 2*dz
	dzs = dz*dz
	dx2 = 2*dx
	dxs = dx*dx
	dxdt = dx*dz
	dxdt4 = dxdt*4

	n2c = 0*ones(nn-2,mm-2)
	n2zc = 0*ones(nn-2,mm-2)
	m2c = 0*ones(nn-2,mm-2)
	m2xc = 0*ones(nn-2,mm-2)
	m2zc = 0*ones(nn-2,mm-2)

	for m = 2:mm-3
       		n2c[:,m] = (rho[:,m+1]-rho[:, m-1])/dz2
       		n2zc[:,m] = (rho[:,m+1] - 2*rho[:,m] + rho[:,m-1])/dzs
	end

    	n2c[:,mm-2] = (rho[:,mm-2]-rho[:,mm-3])/dz
    	n2zc[:,mm-2] = 0*rho[:,mm-2]
    
    	n2c[:,1] = (rho[:,2]-rho[:,1])/dz
    	n2zc[:,1] = 0*rho[:,1]
    
    	n2c = n2c*ccc
    	n2zc = n2zc*ccc
    
	for n = 2:nn-3
    		m2c[n,:] = (rho[n+1,:]-rho[n-1,:])/dx2
    		m2xc[n,:] = (rho[n+1,:] - 2*rho[n,:] + rho[n-1,:])/dxs
    		for m = 2:mm-3
        		m2zc[n,m] = (rho[n+1,m+1] - rho[n+1,m-1] - rho[n-1,m+1] + rho[n-1,m-1])/dxdt4
    		end
    		m2zc[n,1] = 2*(rho[n+1,2] - rho[n+1,1] - rho[n-1,2] + rho[n-1,1])/dxdt4
    		m2zc[n,mm-2] = 2*(rho[n+1,mm-2] - rho[n+1,mm-3] - rho[n-1,mm-2] + rho[n-1,mm-3])/dxdt4
	end
  
    	m2c[1,:] = (rho[2,:] - rho[1,:])/dx
    	m2xc[1,:] = 0*rho[1,:]
    	m2zc[1,2:mm-3] = 2*(rho[2,3:mm-2] - rho[2,1:mm-4] - rho[1,3:mm-2] + rho[1,1:mm-4])/dxdt4
    	m2zc[1,1] = (rho[2,2] - rho[1,2] - rho[2,1] + rho[1,1])/dxdt
    	m2zc[1,mm-2] = (rho[2,mm-2] - rho[2,mm-3] - rho[1,mm-2] + rho[1,mm-3])/dxdt
	
   	m2c[nn-2,:] = (rho[nn-2,:] - rho[nn-3,:])/dx
    	m2xc[nn-2,:] = 0*rho[nn-2,:]
    	m2zc[nn-2,2:mm-3] = 2*(rho[nn-2,3:mm-2] - rho[nn-2,1:mm-4] - rho[nn-3,3:mm-2] + rho[nn-3,1:mm-4])/dxdt4
    	m2zc[nn-2,1] = (rho[nn-2,2] - rho[nn-3,2] - rho[nn-2,1] + rho[nn-3,1])/dxdt
    	m2zc[nn-2,mm-2] = (rho[nn-2,mm-2] - rho[nn-2,mm-3] - rho[nn-3,mm-2] + rho[nn-3,mm-3])/dxdt

    	m2c = ccc*m2c
    	m2xc = ccc*m2xc
    	m2zc = ccc*m2zc
    
# convert to sigma coordinates

	m2 = jHWavesZtoSig(m2c)
	m2x = jHWavesZtoSig(m2xc)
	m2z = jHWavesZtoSig(m2zc)
	n2 = jHWavesZtoSig(n2c)
	n2z = jHWavesZtoSig(n2zc)


	fs = f*ones(size(vmx)) + vmx - (m2.*m2)./(f*n2)
	fstest = minimum(minimum(fs/f))

	ii = findall(n2 -> n2 < 0, n2)


#  	     Make plots of density, mean v

	xgr = 0*ones(nn-2,mm-2)
	zgr = 0*xgr
	zzgr = 0*zgr
	for n = 1:nn-2
    		xtemp = dx*(n-1)
    		xgr[n,:] = xtemp*ones(1,mm-2)/1e5
    		zgr[n,:] = h[n]*(-1*ones(mm-2) + dt*(0:mm-3))/100
		zzgr[n,:] = hmax*(-1*ones(mm-2) + dt*(0:mm-3))/100
	end

	zd = maximum(h)
	z = zd*(-1*ones(mm-2) + dt*(0:mm-3))/100
	x = dx*(0:nn-3)*1e-5
	hh = h/100

	nqq = length(h)
	xff = dx*(0:(nqq-1))*1e-5
	ybb = -maximum(hh)*ones(nqq)
	xffm = maximum(xff)

	fig2 = Figure()
	ax1 = Axis(fig2[1,1])	
	ax2 = Axis(fig2[2,1])
	

	contour!(ax1, x,z,1000*rho, levels = 10, color = :black)
	
	if isempty( ii) == false
		scatter!(ax1,vec(xgr[ii]),vec(zgr[ii]), color = :red, marker = :cross)
	end
	if (fstest/f) < 0
	#		fs changes sign somewhere!
		contour!(ax1,xgr,zgr,fs*1e5,color = :red, levels = 5)
		println(' ')
		println("fs/f becomes negative somewhere: symmetric instability is possible")
		println("	f* X 10^5 is contoured in red along with density")
		println(' ')
		ax1.title = "Density (sigma-t units) and f* X10^5 (1/sec)"
	else
		ax1.title = "Density (sigma-t units)"
	end
	band!(ax1,xff,ybb,-hh,color = :blue)
	ax1.xlabel = "x [km]"
	ax1.ylabel = "z [m]"
	xlims!(ax1,0, xmax)
	ylims!(ax1,-zd/100, 0)
	hidedecorations!(ax1, grid = true, ticklabels = false, label = false)




	band!(ax2,xff,ybb,-hh,color = :blue)
	if abs(vzero) < 0.001
		text!(ax2,0.6*xffm,-zd*0.6/100, text = "v_0 = 0")
	else
		contour!(ax2,xgr,zzgr,vtemp, color = :black, levels = 5)
	end
	ax2.xlabel = "x [km]"
	ax2.ylabel = "z [m]"
	tstr = @sprintf("Mean Alongshore velocity (vzero = %.1f cm/sec)", vzero)
	ax2.title = tstr
	xlims!(ax2,0, xmax)
	ylims!(ax2,-zd/100, 0)

	band!(ax2,xff,ybb,-hh,color = :blue)
	hidedecorations!(ax2, grid = true, ticklabels = false, label = false)

	
	display(GLMakie.Screen(),fig2)


	if ipause > 0.5
		println(' ')
		println("Sleeping for 10 seconds")
		println(' ')
		sleep(10)
	else
		sleep(0.01)
	end
	

	iflag = 0
	if isempty(ii) == false
    		iflag = 1 
	end

	rtup = (rho, iflag)
	return rtup
end								# end of function






"""
	jHWavesZtoSig

	Convert a field from z to sigma coordinates


K. Brink 11/25/2025
"""



function jHWavesZtoSig(gridz)

#   convert gridz (on x,z grid) to gridsig (on x,sigma grid)
#   works for arrays like vm that are (nn-2,mm-2) 
#   Call as
#		gridisg = jHWavesZtoSig(gridz)
#
#  K.H. Brink 11/25/2025  based on Matlab of 5/24/2019

	global nn, mm, dt, h

	th = -1*ones(mm-2)  + dt*(0:(mm-3))	
	zin = maximum(h)*th
	gridsig = NaN*ones(nn-2,mm-2)

	for n = 1:(nn-2)
    		z = th*h[n]
    		zzz = zin
    		if z[end] > zzz[end]
       		 	zzz[end] = z[end]
    		end
		itp = interpolate((vec(zzz),),vec(gridz[n,:]), Gridded(Linear()))
		gridsig[n,:] = itp(z)
	end
	return gridsig
end							# end of function






"""
	jHWavewsSigtoZ(gridsig)

	Convert sigma to z coordinates

K.H. Brink 12/9/2025
"""



function jHWavesSigtoZ(gridsig)

#   Convert gridded field to z coordinates from sigma coordinates
#	call as
#		gridz = jHWavesSigtoZ(gridsig)
#
#	Not used, but provided for completeness
#	K.H. Brink 12/9/2025, based on Matlab of 5/5/03

	global nn, mm, dt, dx, h

	th = -1*ones(mm-2)  + dt*(0:mm-3)
	zout = maximum(h)*th

	gridz = NaN*ones(nn-2,mm-2)

	for n = 1:nn-2
    		z = h[n]*th
    		jj = findall(zout -> zout >= -h(n) && zout <= 0, zout)
		itp = interpolate((vec(z),),vec(gridsig[n,:]),Gridded(Linear()))
		gridz[n,jj] = itp(zout[jj])

    		nss = maximum(size(jj))
    		if nss <(mm-2)
        		jj = find((-1.2*h[n]) < zout &  zout < -h[n]|(zout < -h[n] && zout > -0.2*maximum(h)))
        		gridz[n,jj] = gridsig[n,1]*ones(size(gridz[n,jj]))
    		end  
	end
	return gridz
end								# end of function









"""
	jHWavesCDrvr(ww,yyyy)

	Assemble the matrix to be solved for pressure

K.H. Brink 11/28/2025
"""

function  jHWavesCDrvr(ww,yyyy)

#   Create the big array that gets solved for pressure
#   BB is the pressure field in response to the arbitrary forcing at frequency ww
#	call as
#		rrr = jHWavesCDrvr(ww,yyyy)
#	where yyyy is a dummy
#
#   K. H. Brink 11/28/2025 based on Matlab of 6/27/2018

	global nn, mm, f, rl, dx, dt, h, BB, n2

	w = ww[1] + im*ww[2]
	f2 = f*f
	wf2 = w*f2

	nm = nn*mm
	nrank = nn*mm


#   Do the corners: multiple (commented out) choices provided for the user's taste

#       1,1 (lower, onshore): several options are possible. Uses pzz = 0
     	
	aa1 = 1.
	#   pz = 0
		# aa2 = 0
    		# aa3 = -1
	#   pzz = 0
 		aa2 = -2.
		aa3 = 1.
	#   pxx = 0
    		#AA(1,2*mm+1) = 1
    		#AA(1,mm+1) = -2
	#   pxz = 0
    		#AA(1,3) = -1
    		#AA(1,2*mm+1) = -1
    		#AA(1,2*mm+3) = 1
    		#
    		#AA(1,mm+1) = -1
    		#AA(1,mm+2) = 1
    		#AA(1,2) = -1


    	#
	aa4 = [aa1,aa2,aa3]*wf2/n2[1,1]
	AA = sparse([1,1,1],[1,2,3],aa4,nm,nm)
    	#
	
	

#   upper, onshore
    	#AA[mm,mm] = 1
	aa1 = 1.
	#   pz = 0
    		aa2 = 0
		aa3 = -1.
	#   pzz = 0
    		#AA[mm,mm-1] = -2
    		#AA[mm,mm-2] = 1
	#   pxx = 0;
    		#AA[mm,2*mm] = -2
    		#AA[mm,3*mm] = 1
	#   pxz = 0
    		#AA[mm,mm-2] = -1
    		#AA[mm,3*mm] = -1
    		#AA[mm,3*mm-2] = 1     
    	#
    	
	aa4 = [aa3, aa2, aa1]*wf2/n2[1,mm-2]
	AA = AA + sparse([mm,mm,mm],[(mm-2), (mm-1),(mm)],aa4,nm,nm)
    	

#       Offshore, bottom
    	mr = (nn-1)*mm +1
    	#AA[mr,mr] = 1
	aa1 = 1.
	#   pz = 0
		#aa2 = 0.
    		#aa3 = -1.
	#   pzz = 0
		aa2 =-2.
		aa3 = 1.

	#   pxx = 0
   		#AA[mr,mr-2*mm] = 1
    		#AA[mr,mr-mm] = -2
	#   pxz = 0;
    		#AA[mr,mr-2*mm] = -1
    		#AA[mr,mr + 2] = -1
    		#AA[mr,mr-2*mm + 2] = 1
   	#

	
	aa4 = [aa1, aa2, aa3]*wf2/n2[nn-2,1]
	AA = AA + sparse([mr, mr, mr],[mr, mr+1, mr+2],aa4,nm,nm)
    	#

#       Offshore, top
    	#AA[nm,nm] = 1
	aa1 = 1.

	#   pz = 0
		aa2 = 0.
		aa3 = -1.
	#   pzz = 0
    		#AA(nm,nm-1) = -2
    		#AA(nm,nm-2) = 1
	#   pxx = 0
    		#AA(nm,nm-2*mm) = 1
    		#AA(nm,nm-mm) = -2
	#   pxz = 0
    		#AA(nm,nm-2*mm) = -1
    		#A(nm,nm-2) = -1;
    		#AA(nm,nm-2*mm - 2) = 1
    	#
    	
	aa4 = [aa3, aa2, aa1]*wf2/n2[nn-2,mm-2]
	AA = AA + sparse([nm, nm, nm],[nm-2, nm-1, nm],aa4,nm,nm)
    	#
    
   
#   Do the coastal boundary

	cc = jHWavesCCc(w)


	jl = 2
	jh = mm-1
	ipts = jl:jh							# row

	vv = vec(cc[:,2])
	AA = AA + sparse(vec(ipts),vec(ipts),vv,nm,nm)

	id = ipts -1*ones(mm-2)						# column
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,1]),nm,nm)

	id = ipts + 1*ones(mm-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,3]),nm,nm)

	id = ipts + (mm -1)*ones(mm-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,4]),nm,nm)

	id = ipts + mm*ones(mm-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,5]),nm,nm)

	id = ipts + (mm + 1)*ones(mm-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,6]),nm,nm)

	id = ipts + (2*mm - 1)*ones(mm-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,7]),nm,nm)

	id = ipts + 2*mm*ones(mm-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,8]),nm,nm)

	id = ipts + (2*mm + 1)*ones(mm-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,9]),nm,nm)


  	ibot = (mm+2)*ones(mm-2)  
	AA = AA + sparse(vec(ipts),vec(ibot),vec(cc[:,12]),nm,nm)

	id = ibot - 1*ones(mm-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,11]),nm,nm)

	id = ibot + 1*ones(mm-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,13]),nm,nm)

	id = ibot -mm*ones(mm-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,10]),nm,nm)

	id = ibot + mm*ones(mm-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,14]),nm,nm)

	cc = nothing



#   Do the offshore boundary
#   This assumes the bottom is flat at x = xmax if it is open


	cc = jHWavesCOff(w)
	nlow = nm - mm + 2
	nhigh = nm -1
	ipts =  nlow:nhigh

	AA = AA + sparse(vec(ipts),vec(ipts),vec(cc[:,8]),nm,nm)

	id = ipts -(2*mm +1)*ones(mm-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,1]),nm,nm)

	id = ipts - 2*mm*ones(mm-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,2]),nm,nm)

	id = ipts -(2*mm - 1)*ones(mm-2)
	AA = AA + sparse(vec(ipts),vec(id), vec(cc[:,3]),nm,nm)

	id = ipts -(mm +1)*ones(mm-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,4]),nm,nm)

	id = ipts -mm*ones(mm-2)
	AA = AA + sparse(vec(ipts), vec(id), vec(cc[:,5]),nm,nm)

	id = ipts -(mm - 1)*ones(mm-2)	
	AA = AA + sparse(vec(ipts), vec(id), vec(cc[:,6]),nm,nm)

	id = ipts - 1*ones(mm-2)
	AA = AA + sparse(vec(ipts), vec(id), vec(cc[:,7]),nm,nm)

	id = ipts + 1*ones(mm-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,9]),nm,nm)


	ibot = (nm-2*mm + 2)*ones(mm-2)

	AA = AA + sparse(vec(ipts),vec(ibot),vec(cc[:,12]),nm,nm)
	
	id = ibot -1*ones(mm-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,11]),nm,nm)
	
	id = ibot + 1*ones(mm-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,13]),nm,nm)
	
	id = ibot -mm*ones(mm-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,10]),nm,nm)

	id = ibot + mm*ones(mm-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,14]),nm,nm)
	
	cc = nothing

#   Do the bottom boundary condition

	cc = jHWavesCBot(w)

	ipts = (mm+1):mm:(nrank-2*mm+1)

	AA = AA + sparse(vec(ipts),vec(ipts),vec(cc[:,4]),nm,nm)

	id = ipts - mm*ones(nn-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,1]),nm,nm)

	id = ipts -(mm - 1)*ones(nn-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,2]),nm,nm)

	id = ipts -(mm - 2)*ones(nn-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,3]),nm,nm)

	id = ipts + 1*ones(nn-2)
	AA = AA + sparse(vec(ipts),vec(id), vec(cc[:,5]),nm,nm)

	id = ipts + 2*ones(nn-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,6]),nm,nm)

	id = ipts + mm*ones(nn-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,7]),nm,nm)

	id = ipts + (mm + 1)*ones(nn-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,8]),nm,nm)

	id = ipts + (mm + 2)*ones(nn-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,9]),nm,nm)


	cc = nothing

#   Do the surface boundary condition

	cc = jHWavesSbc(w)

	ipts = (2*mm):mm:(nrank-mm)

	AA = AA + sparse(vec(ipts),vec(ipts) , vec(cc[:,4]),nm,nm)

	id = ipts - (mm + 1)*ones(nn-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,1]),nm,nm)

	id = ipts - 2*ones(nn-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,2]),nm,nm)

	id = ipts - 1*ones(nn-2)
	AA = AA + sparse(vec(ipts),vec(id), vec(cc[:,3]),nm,nm)

	id = ipts + (mm - 1)*ones(nn-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,5]),nm,nm)

	cc = nothing


#   Do the interior points

	for n = 2:(nn-1)

    		cc = jHWavesMatco(w,n)
    		ipts = ((n-1)*mm + 2):1:((n*mm)-1)


		AA = AA + sparse(vec(ipts),vec(ipts),vec(cc[:,5]),nm,nm)

		id = ipts - (mm+1)*ones(mm-2)
		AA = AA + sparse(vec(ipts), vec(id),vec(cc[:,1]),nm,nm)

		id = ipts -mm*ones(mm-2)
		AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,2]),nm,nm)

		id = ipts -(mm-1)*ones(mm-2)
		AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,3]),nm,nm)

		id = ipts -1*ones(mm-2)
		AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,4]),nm,nm)

		id = ipts + 1*ones(mm-2)
		AA = AA + sparse(vec(ipts),vec(id), vec(cc[:,6]),nm,nm)

		id = ipts + (mm - 1)*ones(mm-2)
		AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,7]),nm,nm)

		id = ipts + mm*ones(mm-2)
		AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,8]),nm,nm)

		id = ipts + (mm + 1)*ones(mm-2)
    		AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,9]),nm,nm)
    
		cc = nothing
	end
 

	B = 0*ones(nm)

	inn = Int(trunc(nn/2))
	iforce = mm*inn:mm:(mm*inn + 4*mm)
	B[iforce] = ones(size(iforce))

	bb = AA\B

	area = dx*(nn-3)*maximum(h)
	rrr = jHWavesInt(bb)
	BB = bb*sqrt(rrr*area)
	bref = BB[2*mm-1]
	fact = conj(bref)/abs(bref)
	BB = BB*fact

	@printf("wr, rrr = %.4e %.4eim  %.4e \n", real(w), imag(w), rrr)
	

	return rrr
end							# end of function



"""
	jHWavesCCc

	Compute coefficients for coastal boundary condition
		for use with jHWavesC

K. Brink
"""


function jHWavesCCc(w)

#   Compute coefficients for coastal boundary condition (avoiding very top and bottom grid points)
#		For use with jHWavesC			
#   Modified version to do the closed boundary condition more accurately
#
#	K. Brink 3/10/2026, based on Matlab of 6/11/2018

	global rl, m2, n2, f, vm, h, r, nn, mm, vmx, vmxx, icbc, dx, dt, hx, m2x, m2z

	rr = r[1]
	hh = h[1]

	th = -1*ones(mm-2) + dt*(0:(mm-3))

	wp = w*ones(mm-2) + rl*vm[1,:]
	wps = wp.*wp
	m2n2 = m2[1,:]./n2[1,:]
	fs = (f*ones(mm-2) + vmx[1,:]) - m2[1,:].*m2n2/f
	fsw = f*fs - wps
	fswx = f*vmxx[1,:] - 2*m2x[1,:].*m2n2 + (m2z[1,:].*m2n2).*m2n2 - 2*rl*wp.*vmx[1,:]
	m2n2x = m2x[1,:]./n2[1,:] - (m2z[1,:].*m2n2)./n2[1,:]


	thx = -(hx[1]/hh)*th
	thz = 1/hh

	fsb = fs[1]
	wpsb = wps[1]
	m2n2b = m2n2[1]
	wpb = wp[1]
	fswb = fsw[1]
	fmul = fsw./fswb

	cc = 0*ones(mm-2,14)*(1 + im*1)

	dxdt = dx/dt
	dtdx2 = 2*dt*dx
	dt2 = 2*dt
	dxs = dx*dx


	if icbc < 0.5
# 		      Closed boundary at x = 0

    		d1 = -im*wp*hh
    		d2 =   im*wp.*m2n2*hh
    		d3 = -im*rl*f*hh*ones(mm-2,1)

    		d4 =  -(rr/fswb)*(f*fsb + wpsb)*fmul
    		d5 = (rr/fswb)*(2*wpsb*m2n2b)*fmul
    		d6 =  -(rr/fswb)*(2*wpb*rl*f)*fmul
     

    		C1 = d1
    		C2 = (d1.*thx + d2*thz)
    		C3 = d3
    
    		D4 = d4
	    	D5 = d4*thx[1] + d5*thz
    		D6 = d6
    
    
    		cc[:,2] = -im*C1
    		cc[:,4] = -im*C2*dxdt
    		cc[:,5] = im*C3*2*dx
    		cc[:,6] = im*C2*dxdt
    		cc[:,8] = im*C1
    
    		cc[:,10] = -im*D4
    		cc[:,11] = -im*D5*dxdt
    		cc[:,12] = im*2*dx*D6
    		cc[:,13] = im*D5*dxdt
    		cc[:,14] = im*D4
    
    
    
    #   Now scale it for consistency

     		cc  = cc*hh/(dxdt*dxdt)


 	else 
#		       Open boundary at x = 0          Requires a flat bottom at x = 0

    		b1 = -wp
    		b2 = thz*wp.*m2n2/h[1]
    		b3 = -rl*(f*ones(mm-2) + vmx[1,:]) + (wp.*fswx)./fsw
    		b4 = thz*(rl*vmx[1,:].*m2n2 - ((wp.*m2n2).*fswx)./fsw + wp.*m2n2x)
    		b5 = rl*f*fswx./fsw
    
    		cc[:,1] = 0.5*b2/dtdx2
    		cc[:,2] = b1/dxs - b3*0.5/dx
    		cc[:,3] = -cc[:,1]
    		cc[:,4] = -b4/dt2
    		cc[:,5] = -2*b1/dxs + b5
    		cc[:,6] = b4/dt2
   		cc[:,7] = -cc[:,1]
    		cc[:,8] = b1/dxs + b3*0.5/dx
    		cc[:,9] = cc[:,1]
    
    #   Scale for consistency
    		cc = cc*hh*hh*dt*dt
   
	end								
	return cc
end								# end of function


     


							


"""
	jHWavesCOff(w)

	Compute coefficient for the offshore boundary condition

K. Brink 3/11/2026
"""

function jHWavesCOff(w)

#   Compute the coefficients for the offshore boundary condition
#	for use with jHWavesC
#
#  K. Brink 3/11/2026, baqsed on Matlab of 6/11/2018

	global rl, m2, n2, f, vm, h, hx, nn, mm, vmx, vmxx, m2x, m2z, dx, dt, iobc, r


	hh = h[nn-2]	
	rr = r[nn-2]

	th = -1*ones(mm-2) + dt*(0:(mm-3))

	wp = w*ones(mm-2) + rl*vm[nn-2,:]
	wps = wp.*wp
	m2n2 = m2[nn-2,:]./n2[nn-2,:]
	fs = (f*ones(mm-2) + vmx[nn-2,:]) - m2[nn-2,:].*m2n2/f
	fsw = f*fs - wps
	mmmm = (m2z[nn-2,:].*m2n2).*m2n2
	fswx = f*vmxx[nn-2,:] - 2*rl*wp.*vmx[nn-2,:] -2*m2x[nn-2,:].*m2n2 + mmmm
	m2n2x = m2x[nn-2,:]./n2[nn-2,:] - (m2n2.*m2z[nn-2,:])./n2[nn-2,:]


	thz = 1/hh
	cc = 0*ones(mm-2,14)*(1 + im*1)

	if iobc > 0.5
#  			 du/dx = 0 (Open) boundary condition offshore
    
    
		b1 = -wp
		b2 = thz*wp.*m2n2
		b3 = -rl*( f*ones(mm-2) + vmx[nn-2,:]) + wp.*(fswx./fsw)
		b4 = thz*(rl*vmx[nn-2,:].*m2n2 - ((wp.*m2n2).*fswx)./fsw  + wp.*m2n2x)
		b5 = rl*f*fswx./fsw

		dxdt2 = 4*dx*dt
		dx2 = 2*dx
		dxs = dx*dx



		cc[:,1] = b2/dxdt2
		cc[:,2] =  b1/dxs - b3/dx2
		cc[:,3] = -cc[:,1]
		cc[:,4] = -b4*0.5/dt
		cc[:,5] = b5 - 2*b1/dxs
		cc[:,6] = -cc[:,4]
		cc[:,7] = -cc[:,1]
		cc[:,8] = b1/dxs + b3/dx2
		cc[:,9] = cc[:,1]

		cc = cc*hh*hh*dt*dt

	else
    # 			  No flow through the offshore boundary
    
    
		thx = -(hx[nn-2]/hh)*th
		dxdt = dx/dt

		fsb = fs[1]
		wpsb = wps[1]
		m2n2b = m2n2[1]
		wpb = wp[1]
		fswb = fsw[1]
		fmul = fsw./fswb


   		 d1 = -im*wp*hh
   		 d2 =  im*wp.*m2n2*hh
 		 d3 = -im*rl*f*hh*ones(mm-2,1)

  		 d4 =  -(rr/fswb)*(f*fsb + wpsb)*fmul
 		 d5 = (rr/fswb)*(2*wpsb*m2n2b)*fmul
 		 d6 =  -(rr/fswb)*(2*wpb*rl*f)*fmul
     

 		 C1 = d1
 		 C2 = (d1.*thx + d2*thz)
    		 C3 = d3
    
   		 D4 = d4
    		 D5 = d4*thx[1] + d5*thz
    		 D6 = d6
    
    
    		 cc[:,2] = -im*C1
    		 cc[:,4] = -im*C2*dxdt
    		 cc[:,5] = im*C3*2*dx
    		 cc[:,6] = im*C2*dxdt
   		 cc[:,8] = im*C1
    
    		 cc[:,10] = -im*D4
    		 cc[:,11] = -im*D5*dxdt
    		 cc[:,12] = im*2*dx*D6
    		 cc[:,13] = im*D5*dxdt
    		 cc[:,14] = im*D4
    
    
    
    #   Now scale it for consistency

     		cc  = cc*hh/(dxdt*dxdt)
 	end
	return cc
end								# end of function




"""
	jHWavesCBot

	Compute coefficient for the bottom boundary condition
		for use with jHWavesC

K. Brink	3/11/2026
"""	




function  jHWavesCBot(w)
#   Compute coefficients for bottom boundary condition
#		for use with jHWavesC
#
# K. Brink 3/11/2026, based on Matlab of 6/11/2018

	global nn, mm, dx, dt, r, rx, m2, m2x, m2z, n2, n2z, f, rl, vm, vmx, vmxx, h, hx, hxx


	m2n2 = (m2[:,1]./n2[:,1])
	fs = f*ones(nn-2) + vmx[:,1] - m2[:,1].*m2n2/f
	fsx = vmxx[:,1] -2*m2x[:,1].*m2n2/f + (m2n2.*m2z[:,1]).*m2n2/f
	wp = w*ones(nn-2) + rl*vm[:,1]
	wps = wp.*wp
	fsw = f*fs - wps
	fspw = f*fs + wps
	rfsw = r./fsw
	fswx = f*fsx -2*rl*vmx[:,1].*wp
	m2n2x = m2x[:,1]./n2[:,1] - (m2n2.*m2z[:,1])./n2[:,1]
	fp = f*ones(nn-2) + vmx[:,1]


	th = -1
	thz = ones(nn-2,1)./h
	thxz = -(hx./(h.*h))
	thx = (hx./h)
	thxbx = (hxx./h) - thx.*thx

	m2n2bx = m2n2x -hx.*(m2z[:,1]./n2[:,1] - (n2z[:,1].*m2n2)./n2[:,1])
	fsbx = vmxx[:,1] - (hx.*m2x[:,1]/f) - (m2n2bx.*m2[:,1]/f) - (m2n2/f).*(m2x[:,1] - hx.*m2z[:,1])
	vmbx = vmx[:,1] - hx.*m2[:,1]/f
	fswbx = f*fsbx - 2*rl*wp.*vmbx

	rrbx = (rx./fsw)  - 2*(r.*fswbx)./(fsw.*fsw)
	coef = 2*rl*(wp.*vmbx).*m2n2 + wps.*m2n2bx

	b1 = -rfsw.*fspw
	b2 = rfsw.*(2*wps.*m2n2)
	b3 = im*wp.*m2n2 -im*wp.*hx - fspw.*rrbx + rfsw.*(-2*rl*vmbx.*wp - f*fsbx + 2*rl*fs.*wp)
	b4 = -2*f*rl*(wp.*rfsw)
	b5 = -im*(wp./n2[:,1]).*(f*fp - wp.*wp)  + im*(wp.*m2n2).*hx
	b5 = b5 + 2*(wps.*m2n2).*rrbx + rfsw.*(-(rl/f)*(wp.*m2n2).*fspw + 2*coef)
	b6 = im*rl*f*m2n2 - im*rl*f*hx - 2*rl*f*wp.*rrbx + rl*rl*rfsw.*(-2*f*vmbx + fspw)


	d1 = b1
	d2 = b1.*thx + b2.*thz
	d3 = b1.*thxbx + b3.*thx + b5.*thz + b2.*thxz
	d4 = b3 + b4
	d5 = b6


#   Scale for consistency
	d1 = dt*d1./(thz)
	d2 = dt*d2./(thz)
	d3 = dt*d3./(thz)
	d4 = dt*d4./(thz)
	d5 = dt*d5./(thz)

	dxdt4 = 4*dx*dt
	dxs = dx*dx
	dx2 = 2*dx
	dt2 = 2*dt

	cc = 0*ones(nn-2,9)*(1 + im*1)
	cc[:,1] = d2/dxdt4
	cc[:,2] = d1/dxs - d4/dx2
	cc[:,3] = -d2/dxdt4
	cc[:,4] = -d3/dt2
	cc[:,5] = -2*d1/dxs + d5
	cc[:,6] = d3/dt2
	cc[:,7] = -d2/dxdt4
	cc[:,8] = d1/dxs + d4/dx2
	cc[:,9] = d2/dxdt4

	cc = im*cc

	return cc
end								# end of function



"""
	jHWavesSbc(w)

	Compute coefficients for surface boundary condition

K.H. Brink  11/28/2025
"""


function jHWavesSbc(w)

#   Compute coefficients for the surface boundary condition
#	call as
#		cc = jHWavesSbc(w)
#
#  K.H. Brink 11/28/2025 based on Matlab of 12/23/2004


	global del, nn, mm, f, rl, vm, m2, n2, vmx, h, dx, dt

	g = 980.

	cc = 0*ones(nn-2,5)*(1 + im*1)

	wp = w*ones(nn-2) + rl*vm[:,mm-2]
	m2n2 = m2[:,mm-2]./n2[:,mm-2]
	fs = f*ones(nn-2) + vmx[:,mm-2] - m2[:,mm-2].*m2n2/f
	fsw = f*fs - wp.*wp

	thz = ones(nn-2)./h
	thx = 0.0

	s1 = ((-wp./n2[:,mm-2]).*(fsw + m2[:,mm-2].*m2n2))./thz
	s2 = (wp.*m2n2)./thz
	s3 = (rl*f*m2n2 - (wp*del).*fsw/g)./thz

	c1 = s1.*thz + s2*thx

	cc[:,1] = -s2
	cc[:,2] = -c1*dx/dt
	cc[:,3] = 2*dx*s3
	cc[:,4] = c1*dx/dt
	cc[:,5] = s2

	cc = cc*dt/dx
	return cc
end						# end of function





"""
	jHWavesMatco

	Compute coefficients for the model interior

K.H. Brink 11/29/2025
"""

 
function jHWavesMatco(w,nq)

#   Get coefficients for the big array to be solved 
#   This works for a particular offshore grid point n
#	call as
#		cc = jHWavesMatco(w,nq)
#
#	K.H. Brink 11/29/2025


	global nn, mm, dt, dx, f, rl, vm, vmx, vmxx, m2, n2, n2z, m2z, m2x, h, hx, hxx


	n = nq-1

	wp = w*ones(mm-2)+ rl*vm[n,:]

	wps = wp.*wp
	m2n2 = (m2[n,:]./n2[n,:])

	fs = f*ones(mm-2) + vmx[n,:] - (m2n2.*m2[n,:])/f
	fsw = f*fs - wps


	fsx = vmxx[n,:] -2*m2x[n,:].*m2n2/f + ((m2z[n,:].*m2n2).*m2n2)/f
	fswx = f*fsx -2*rl*wp.*vmx[n,:]


	fswz = m2x[n,:]- 2*rl*wp.*m2[n,:]/f - 2*m2z[n,:].*m2n2 + ((m2n2.*m2n2).*n2z[n,:])
	fp = f*ones(mm-2) + vmx[n,:]
	m2n2z = m2z[n,:]./n2[n,:] - (n2z[n,:].*m2n2)./n2[n,:]
	rq = (fswx - m2n2.*fswz)./fsw
	m2n2x = m2x[n,:]./n2[n,:] - (m2z[n,:].*m2n2)./n2[n,:]

	a1 = -wp
	a2 = 2*m2n2.*wp
	a3 = -(wp.*(f*fp - wps))./n2[n,:]
	a4 = wp.*rq + wp.*m2n2z
	a5 = -(wp.*m2n2).*rq + rl*vmx[n,:].*m2n2 + wp.*m2n2x + (rl/f)*m2n2.*(f*f*ones(mm-2) - wps)
	a5 = a5 - (rl*(m2[n,:]/f)./n2[n,:]  - (wp.*n2z[n,:])./(n2[n,:].*n2[n,:])).*(f*fp - wps)
	a5 = a5 - (wp./n2[n,:]).*(m2z[n,:].*m2n2 + m2[n,:].*m2n2z)
	a6 = rl*f*rq + wp*rl*rl + rl*f*m2n2z

	th = -1.0*ones(mm-2) + dt*(0:(mm-3))						
	thz = 1/h[n]
	thx = -hx[n]*th*thz
	thxx = -(hxx[n]/h[n])*th - 2*(hx[n]/h[n])*thx
	thxz = -hx[n]*thz*thz

	b1 = a1
	b2 = 2*thx.*a1 + a2.*thz
	b3 = (thx.*thx).*a1 + a3*thz*thz
	b3 = b3 + thz*thx.*a2
	b4 = a1.*thxx + a2.*thxz + a4.*thx + a5.*thz
	b5 = a4
	b6 = a6

	dt2 = 2*dt
	dxdts2 = 2*dx/(dt*dt)
	dxh = 2/dx
	dxdt = dx/dt


	cc = 0*ones(mm-2,9)*(1 + im*1)
	cc[:,1] = b2/dt2
	cc[:,2] = b1*dxh - b5
	cc[:,3] = -b2/dt2
	cc[:,4] = b3*dxdts2 - b4*dxdt
	cc[:,5] = -4*b1/dx + 2*dx*b6 - 2*b3*dxdts2
	cc[:,6] = b3*dxdts2 + b4*dxdt
	cc[:,7] = -b2/dt2
	cc[:,8] = b1*dxh + b5
	cc[:,9] = b2/dt2

#   Scale for consistency
	cc = cc*h[n]*h[n]*dt*dt/dx
	return cc
end							# end of function







"""
	jHWavesInt(bb)

	Compute the inverse norm of the pressure for each iteration

K.H. Brink 12/1/2025
"""

function jHWavesInt(bb)

#   Compute the inverse integral for the given iteration
#   bb comes in the form of an (nn*mm,1) array
#	call as
#		rrr = jHWavesInt(bb) 
#    
#   K.H. Brink 12/1/2025 based on Matlab of 5/22/2018

	global nn, mm, dt, dx, h

	dz = h[1]*dt
	nlow = mm + 2
	nhigh = 2*mm - 1
	rr = 0.25*dz*(abs(bb[nlow])^2 + abs(bb[nhigh])^2)
	rr = rr+ dz*( 0.5*sum(abs.(bb[(nlow + 1):(nhigh - 1)]).^2))

	dz = h[nn-2]*dt
	ntop = (nn-1)*mm - 1
	nbot = (nn-2)*mm + 2
	rt = 0.25*dz*(abs(bb[ntop])^2 + abs(bb[nbot])^2)
	rt = rt + dz*( 0.5*sum(abs.(bb[(nbot+1):(ntop-1)]).^2))

	rrr = 0

	for n = 3:(nn-2)
    		dz = h[n-1]*dt
    		ntop = n*mm - 1
    		nbot = (n-1)*mm + 2
    		rs = 0.5*dz*(abs(bb[nbot])^2 + abs(bb[ntop])^2)
    		rs = rs + dz*(sum(abs.(bb[(nbot+1):(ntop-1)]).^2))
    		rrr = rrr + rs
	end

	rrr = rrr + rt + rr
	rrr = 1/rrr
	return rrr
end								# end of function






"""
	jHWavesCUvwr(ww,BB)

	Compute the wave's u,v, w,rho, given pressure
		for use with jHWavesC

K.H. Brink 3/11/2026
"""


function jHWavesCUvwr(ww,BB)

#   Compute u, v, w and rho once the pressure is known
#		for use with jHWavesC
#   ww is the frequency and BB is the pressure vector
#	call as
#		utup = jHWavesCUvwr(ww,BB)
#	BB is an nm by 1 array
#
#  K.H. Brink 2/5/2026, based on matlab of 5/22/2018, 3/11/2026



	global nn, mm, dx, dt, rl, f, vm, vmx, n2, m2, h, hx

	w = ww

	u = 0*ones(nn-2,mm-2)*(1 +im)
	v = copy(u)
	wvel = copy(u)
	rho = copy(u)
	pp = copy(u)
	ptop = 0*ones(nn-2)*(1 + im)

	tol = 0.01
	g = 980
	rhob = 1.03
	dt2 = 2*dt
	dx2 = 2*dx
	th = -1*ones(mm-2) + dt*(0:(mm-3))

	for n = 1:(nn-2)
    
    		irow = (n*mm + 2):((n+1)*mm -1)
    		wp = w*ones(mm-2) + rl*vm[n,:]
    		m2n2 = m2[n,:]./n2[n,:]
    		fs = f*ones(mm-2) + vmx[n,:] - m2[n,:].*m2n2/f
    		fsw = f*fs - wp.*wp
    		thx = -(hx[n]/h[n])*th
    
    		ptop[n] = BB[((n+1)*mm-1)]
    		p = BB[irow]
		ip = (irow) + (ones(Int,mm-2))
		imin = (irow) - (ones(Int,mm-2))
    		pt = (BB[ip] - BB[imin])/dt2
    		pxp = (BB[(irow+mm*ones(Int,mm-2))] - BB[(irow-mm*ones(Int,mm-2))])/dx2
       
    		px = pxp + thx.*pt
    		pz = pt/h[n]
		vt = fs.*px + rl*wp.*p - (wp.*wp.*m2n2).*pz/f
		v[n,:] = vt./fsw
		
    		wvel[n,:] = -(m2n2./fsw).*(-im*wp.*px - im*rl*f*p)
    		wtemp = ((wp.*pz)./n2[n,:]).*(1.0*ones(mm-2) +(((m2n2.*m2[n,:])./fsw)))
    		wvel[n,:] = wvel[n,:] - im*(wtemp)

    		rho[n,:] = conj((-pz/g)')

		utemp = -im*wp.*px - im*rl*f*p
		utemp = utemp + im*(wp.*m2n2).*pz
		u[n,:] = utemp./fsw

    		pp[n,:] = p
	end


	u = u/rhob
	v = v/rhob
	wvel = wvel/rhob

	utup = (u, v, wvel, rho)


#   Plot results


	xgr = 0*ones(nn-2,mm-2)
	zgr = 0*ones(nn-2,mm-2)
	for n = 1:nn-2
   		xtemp = dx*(n-1)
    		xgr[n,:] = xtemp*ones(1,mm-2)/1e5
    		zgr[n,:] = h[n]*(-1*ones(mm-2) + dt*(0:mm-3))/100
	end

	hh = -h/100
	xpl = dx*(0:nn-3)/1e5
	ax = [0 maximum(xpl) minimum(hh) 0]
	z = th*maximum(h)/100
	xf = [0 xpl' maximum(xpl) 0]
	zf = [minimum(hh) hh' minimum(hh) minimum(hh)]
	xdep = vec(xgr[:,1])
	xmax = maximum(xdep)
	hmax = maximum(h)


	fig5 = Figure()
	

	colsize!(fig5.layout, 1, Aspect(1, 1.5))
				
	ax1 = Axis(fig5[1,1], title = "uplot")				# top panel
	ax3 = Axis(fig5[2,1], title = "wplot")				# lower
	ax2 = Axis(fig5[1,2], title = "vplot")
	ax4 = Axis(fig5[2,2], title = "rplot")


	rmm = maximum(maximum(abs.(real(pp))))
	amm = maximum(maximum(abs.(imag(pp))))
	rata = amm/rmm
	pmult = 1.
	if rata < tol
		pmult = NaN
	end


#		plot u
	rmm = maximum(maximum(abs.(real(u))))
	amm = maximum(maximum(abs.(imag(u))))
	rata = rmm/amm
	umult = 1.
	if rata < tol
		umult = NaN
	end

	civ = jHWavesCCon(u,10)
	
	zmax = maximum(-hh)
	yupper = hh
	nzz = length(yupper)
	ylower = 0*yupper -zmax*ones(nzz)
	band!(ax1,xdep,ylower,yupper,color = (:blue,0.7))
	lines!(ax1,xdep,hh, color = (:blue,0.7))
	xlims!(ax1,0,xmax)
	ylims!(ax1,-hmax/100,0)
	hidedecorations!(ax1, grid = true, ticklabels = false, label = false)
 	ax1.ylabel = "z [m]"
	ax1.xlabel = "x [km]"
	ax1.title = "u"

	contour!(ax1,xgr,zgr,umult*real(u), levels = civ, color = :black)
	contour!(ax1,xgr,zgr,umult*real(u),levels = vec([0 0]), linewidth = 2, color = :black)
	contour!(ax1,xgr,zgr,imag(u), levels = civ, color = :red)
	contour!(ax1,xgr,zgr,imag(u),levels = vec([0 0]), linewidth = 2, color = :red)

	



#		plot v
	civ = jHWavesCCon(v,10)

	
	rmm = maximum(maximum(abs.(real(v))))
	amm = maximum(maximum(abs.(imag(v))))
	rata = amm/rmm
	vmult = 1.
	if rata < tol
		vmult = NaN
	end


	band!(ax2,xdep,ylower,yupper,color = (:blue,0.7))
	lines!(ax2,xdep,hh, color = (:blue,0.7))
	xlims!(ax2,0,xmax)
	ylims!(ax2,-hmax/100,0)
	hidedecorations!(ax2, grid = true, ticklabels = false, label = false)
 	ax2.ylabel = "z [m]"
	ax2.xlabel = "x [km]"
	ax2.title = "v"

	contour!(ax2,xgr,zgr,real(v), levels = civ, color = :black)
	contour!(ax2,xgr,zgr,real(v),levels = vec([0 0]), linewidth = 2, color = :black)
	contour!(ax2,xgr,zgr,vmult*imag(v), levels = civ, color = :red)
	contour!(ax2,xgr,zgr,vmult*imag(v),levels = vec([0 0]), linewidth = 2, color = :red)


#		plot w

	civ = jHWavesCCon(wvel,10)
	

	rmm = maximum(maximum(abs.(real(w))))
	amm = maximum(maximum(abs.(imag(w))))
	rata = rmm/amm
	wmult = 1.
	if rata < tol
		wmult = NaN
	end

	band!(ax3,xdep,ylower,yupper,color = (:blue,0.7))
	lines!(ax3,xdep,hh, color = (:blue,0.7))
	xlims!(ax3,0,xmax)
	ylims!(ax3,-hmax/100,0)
	hidedecorations!(ax3, grid = true, ticklabels = false, label = false)
 	ax3.ylabel = "z [m]"
	ax3.xlabel = "x [km]"
	ax3.title = "w"

	contour!(ax3,xgr,zgr,imag(wvel), levels = civ, color = :red)
	contour!(ax3,xgr,zgr,imag(wvel),levels = vec([0 0]), linewidth = 2, color = :red)
	contour!(ax3,xgr,zgr,wmult*real(wvel), levels = civ, color = :black)
	contour!(ax3,xgr,zgr,wmult*real(wvel),levels = vec([0 0]), linewidth = 2, color = :black)



#		plot rho
	civ = jHWavesCCon(rho,10)

	rmm = maximum(maximum(abs.(real(rho))))
	amm = maximum(maximum(abs.(imag(rho))))
	rata = amm/rmm
	rmult = 1.
	if rata < tol
		rmult = NaN
	end


	band!(ax4,xdep,ylower,yupper,color = (:blue,0.7))
	lines!(ax4,xdep,hh, color = (:blue,0.7))
	xlims!(ax4,0,xmax)
	ylims!(ax4,-hmax/100,0)
	hidedecorations!(ax4, grid = true, ticklabels = false, label = false)
 	ax4.ylabel = "z [m]"
	ax4.xlabel = "x [km]"
	ax4.title = "rho"

	contour!(ax4,xgr,zgr,real(rho), levels = civ, color = :black)
	contour!(ax4,xgr,zgr,real(rho),levels = vec([0 0]), linewidth = 2, color = :black)
	contour!(ax4,xgr,zgr,rmult*imag(rho), levels = civ, color = :red)
	contour!(ax4,xgr,zgr,rmult*imag(rho),levels = vec([0 0]), linewidth = 2, color = :red)



	display(GLMakie.Screen(),fig5)


#       Plot surface values

#		Scale so that abs(v(0,0) = 1
	ssff = maximum(abs.(v[:,end]))
	ratp = 1/ssff;

	fig6 = Figure()
	

	colsize!(fig6.layout, 1, Aspect(1, 1.5))
				
	ax1 = Axis(fig6[1,1], title = "pplot")				# top panel
	ax3 = Axis(fig6[3,1], title = "uplot")				# lower
	ax2 = Axis(fig6[2,1], title = "vplot")
	ax4 = Axis(fig6[4,1], title = "hplot")




#		plot surface p
	lines!(ax1,vec(xdep),real(vec(ptop*ratp)), color = :black)
	lines!(ax1,vec(xdep),pmult*imag(vec(ptop*ratp)), color = :red)

	lines!(ax1,vec(xdep),vec(0*xdep),color = :black)
	xlims!(ax1,0,xmax)
	hidedecorations!(ax1, grid = true, ticklabels = false, label = false)
 	ax1.ylabel = "p [dyne/cm^2]"
	ax1.xlabel = "x [km]"
	ax1.title = "Surface Pressure (all amplitudes arbitrary, but consistent)"


#		plot v
	lines!(ax2,vec(xdep),real(vec(v[:,end]*ratp)), color = :black)
	lines!(ax2,vec(xdep),vmult*imag(vec(v[:,end]*ratp)), color = :red)

	lines!(ax2,vec(xdep),vec(0*xdep), color = :black)
	xlims!(ax2,0,xmax)
	hidedecorations!(ax2, grid = true, ticklabels = false, label = false)
 	ax2.ylabel = "v [cm/s]"
	ax2.xlabel = "x [km]"
	ax2.title = "Surface v"


#		plot u
	lines!(ax3,xdep,umult*real(u[:,end]*ratp), color = :black)
	lines!(ax3,xdep,imag(u[:,end]*ratp), color = :red)

	lines!(ax3,xdep,0*xdep,color = :black)	
	xlims!(ax3,0,xmax)
	hidedecorations!(ax3, grid = true, ticklabels = false, label = false)
 	ax3.ylabel = "u [cm/s]"
	ax3.xlabel = "x [km]"
	ax3.title = "Surface u"

#		plot depth
	
	band!(ax4,xdep,ylower,yupper,color = (:blue,0.7))
	lines!(ax4,xdep,hh, color = (:blue,0.7))
	xlims!(ax4,0,xmax)
	ylims!(ax4,-hmax/100,0)
	
	tt = maximum(maximum(abs.(vm)))
	if tt > 0.1
		contour!(ax4,xgr,zgr,vm, levels = 8, color = :black)
	end

	hidedecorations!(ax4, grid = true, ticklabels = false, label = false)
 	ax4.ylabel = "z [m]"
	ax4.xlabel = "x [km]"
	ax4.title = "Depth profile and mean velocity"
	
	display(GLMakie.Screen(),fig6)

	return utup
end									# end of function





"""
	ci = jHWavesCCon(u,nc)

	define contour intervals, given the range

K. Brink
"""

function jHWavesCCon(u, nc)
	


"""	
	jHWavesCCon(u, nc)

	function to make consistent contour intervals for a complex array u
		for use with jHWavesC
	nc = # of contours
	Returns an array of contour intervals
K. Brink 3/15/2026
"""

	urma = maximum(maximum(real(u)))
	urmi = minimum(minimum(real(u)))
	uima = maximum(maximum(imag(u)))
	uimi = minimum(minimum(imag(u)))
	umi = uimi
	if urmi < uimi
		umi = urmi
	end
	uma = uima
	if urma > uima
		uma = urma
	end
	if umi > 0
		umi = 0
	end	
	if uma < 0
		uma = 0
	end
	ci = (uma-umi)/nc

	cc = umi:ci:uma
	(zz,ind) = findmin(abs.(cc))
	cont = cc .- zz
		
	return cont
end




"""
	jHWavesCEngDiag(war,BB)

	Compute and plot energy diagnosticas for the resulting wave
		for use with jHWavesC

K.H. Brink  3/11/2026
"""


function jHWavesCEngDiag(war,BB)

#  Compute energy diagnistics: most useful for unstable waves
#		for use with jHWavesC
#	call as
#		jHWavesCEngDiag(war,BB)
#
#  K.H. Brink 12/1/2025bade on Matlab of 6/11/2018, 3/11/2026


	global nn, mm, dx, dt, f, rl, del, r, rx, h, hx, hxx, n2, m2, m2x, m2z, n2z, vm, vmx, vmxx

	con1 = 0
	con2 = 0
	con3 = 0
	con4 = 0
	con5 = 0

	eke = 0
	epe = 0
	epes = 0

	g = 980
	rhob = 1.03
	gr = g*rhob

	utup = jHWavesCUvwr(war,BB)
	u = utup[1]
	v = utup[2]
	wvel = utup[3]
	rho = utup[4]

	uc = conj(u)
	vc = conj(v)
	wvelc = conj(wvel)
	rhoc = conj(rho)

	w = war
	m2n2 = m2./n2
	dx2 = 2*dx
	dt2 = 2*dt
	dxs = dx*dx

	urhocov = 0;

	for n = 1:nn-2
    		dz = dt*h[n]
    		dxl = dx
    		if n == 1 
        		dxl = dx/2
    		end
		if n == nn-2
        		dxl = dx/2
		end			

    
    		itop = ((n+1)*mm) - 1
    		epes = epes + dxl*BB[itop]*conj(BB[itop])
    
    		xx = u[n,:].*uc[n,:] + v[n,:].*vc[n,:]
    		yy = (rho[n,:].*rhoc[n,:])./n2[n,:]
    		eke = eke + dxl*dz*(sum(xx) - (xx[1] +xx[mm-2])/2)
    		epe = epe + dxl*dz*(sum(yy) - (yy[1] +yy[mm-2])/2)

    
    		s1 = (vc[n,:].*u[n,:] + uc[n,:].*v[n,:]).*vmx[n,:]
    		s2 = (vc[n,:].*wvel[n,:] + (wvelc[n,:].*v[n,:])).*m2[n,:]/f
    		s3 = wvelc[n,:].*rho[n,:] + wvel[n,:].*rhoc[n,:]
    		s5 = (u[n,:].*rhoc[n,:] + uc[n,:].*rho[n,:]).*m2n2[n,:]
    
    		con1 = con1 + dxl*dz*(sum(s1) - (s1[1] +s1[mm-2])/2)
    		con2 = con2 + dxl*dz*(sum(s2) - (s2[1] +s2[mm-2])/2)
    		con3 = con3 + dxl*dz*(sum(s3) - (s3[1] +s3[mm-2])/2)
    		con5 = con5 + dxl*dz*(sum(s5) - (s5[1] +s5[mm-2])/2)
    
		ssig = m2[n,:]./abs.(m2[n,:])

    		urhocov = urhocov + dxl*dz*sum(u[n,:].*rhoc[n,:].*ssig)
    
   		wp = w + rl*vm[n,1]
    		wps = wp*wp
    		fs = f + vmx[n,1] - m2[n,1]*m2n2[n,1]/f
    		fsw = f*fs - wps
    		rfsw = r[n]/fsw
    
     		imain = n*mm+2					# bottom point
     		px = (BB[imain+mm] - BB[imain-mm])/dx2
   
     		pbx = px
    		 pbxc = conj(pbx)
     
    		 Ue = rfsw*(-im*wp*u[n,1]- f*v[n,1])
     		Uec = conj(Ue)
    
     		con4 = con4 + dxl*(pbx*Uec + pbxc*Ue)
    
 	end



	eke = real(eke*rhob/2)
	epes = real(epes*del/(gr*2))
	epe = real(0.5*epe*g*g/rhob)

	con1 = -real(con1)/2
	con2 = -real(con2)/2
	con3 = -real(con3)*g/2
	con4 = real(con4)/2
	con5 = real(con5)*g/2

	println(' ')
	@printf("EKE = %.4e, Interior EPE = %.4e \n", eke, epe)
	@printf("Surface EPE = %.4e \n", epes)


#   Plot results

	fig10 = Figure()
	

	#colsize!(fig6.layout, 1, Aspect(1, 1.5))
				
	ax1 = Axis(fig10[1,1])	





#				 Define arrows
	arx = [-3, 0, -3]
	ary = [-3, 0, 3]
	arvx = [-2, 0, 2]
	arvy = [-3, 0, -3]

	x = [35, 75, 75, 35, 35]
	y = [60, 60, 90, 90, 60]

	lines!(ax1,x,y, color = :black)
	limits!(ax1,0, 100, 0, 100)

#   text!(ax, 0.5, 0.5, text = "Label", color = :blue, fontsize = 20)

	y = y .- 50
	lines!(ax1,x,y, color = :black)
	text!(ax1, 50, 80, text = "EKE", color = :black, fontsize = 20)
	stt = string(@sprintf("%.3e",eke))
	text!(45,70, text = stt, color = :black,fontsize = 20)
	text!(50,30,text = "EPE", color = :black, fontsize = 20)
	eee = epe + epes
	stt = string(@sprintf("%.3e",eee ))
	text!(45,19, text = stt, color = :black,fontsize = 20)

	x = [15, 35]
	y = [85, 85]
	lines!(ax1,x,y, color = :black)
	if con1 >= 0
 		lines!(ax1,arx .+ 35,ary .+ 85, color = :black)
	else
		lines!(ax1, -arx .+ 15, ary .+ 85, color = :black)
	end
	text!(ax1, 1,90,text = "Barotropic MKE to EKE", color = :black, fontsize = 18)
	stt = string(@sprintf("%.3e",con1))
	text!(ax1,19,80,text = stt, color = :black, fontsize = 15)

	y = y .- 20
	lines!(ax1,x,y, color = :black)
	if con2 >= 0
 		lines!(ax1,arx .+ 35,ary .+ 65, color = :black)
	else
		lines!(ax1, -arx .+ 15, ary .+ 65, color = :black)
	end
	text!(ax1, 5,70,text = "Shear MKE to EKE", color = :black, fontsize = 18)
	stt = string(@sprintf("%.3e",con2))
	text!(ax1,17,60,text = stt, color = :black, fontsize = 15)

	y = [25, 25]
	lines!(ax1,x,y, color = :black)
	if con5 >= 0
 		lines!(ax1,arx .+ 35,ary .+ 25, color = :black)
	else
		lines!(ax1, -arx .+ 15, ary .+ 25, color = :black)
	end
	text!(ax1, 13, 30,text = "MPE to EPE", color = :black, fontsize = 18)
	stt = string(@sprintf("%.3e",con5))
	text!(ax1,19,20,text = stt, color = :black, fontsize = 15)

	x = [55, 55]
	y = [40, 60]
	lines!(ax1,x,y, color = :black)
	if con3 >= 0
 		lines!(ax1,arvx .+ 55,arvy .+ 60, color = :black)
	else
		lines!(ax1, arvx .+ 55, -arvy .+ 40, color = :black)
	end
	text!(ax1, 34,50,text = "EPE to EKE", color = :black, fontsize = 18)
	stt = string(@sprintf("%.3e",con3))
	text!(ax1,58,50,text = stt, color = :black, fontsize = 15)

	x = [75, 95]
	y = [75, 75]
	lines!(ax1,x,y, color = :black)
	if con4 <= 0
 		lines!(ax1,arx .+ 95,ary .+ 75, color = :black)
	else
		lines!(ax1, -arx .+ 75, -ary .+ 75, color = :black)
	end
	text!(ax1, 80,80,text = "Dissipation", color = :black, fontsize = 20)
	stt = string(@sprintf("%.3e",con4))
	text!(ax1,78,70,text = stt, color = :black, fontsize = 15)
	#	This should be a negative number



	stt = string(@sprintf("%.3e %.3eim",real(war), imag(war)))
	ax1.title = string("Energy diagnostics for ω = ", stt, " 1/sec")
		
	#hidedecorations!(ax1, label = false, ticklabels = false, ticks = false )
	hidedecorations!(ax1)

	display(fig10)


end									# end of function




"""
	jHWavesCNorm(pin)

	Normalize the pressure field

K.H. Brink 12/1/2025
"""


function jHWavesCNorm(pin)

#   Normalize pressure mode for complex frequency
#   This is  convenient but not rigorous (does not matter)
#	call as
#		p = jHWavesCNorm(pin)
#	returns a normalized pressure
#
#	K.H. Brink 12/1/2025 based on Matlab of 6/11/2018


	global nn, mm, dx, dt, h, hx, BB, f

	dz = h[1]*dt
	pc = pin[1,:]
	con1 = sum(abs2.(pc)) - 0.5*(abs2(pc[1]) + abs2(pc[mm-2]))
	con1 = con1*dz
	
	pb = hx.*abs2.(pin[:,1])
	con2 = sum(pb) - 0.5*(pb[1] + pb[end])
	con2 = con2*dx

	con = sqrt(abs((con1 + con2)/abs(f)))

	p =  pin/con
	BB = BB/con
	return p

end							# end of function






"""
	jHWavesRtoC(tupR)

	Convert the input tuple for real frequency calculations
		to one for complex frequency
	
K. Brink 3/15/2026
"""

function jHWavesRtoC(tupR)

#	Convert real-frequency input tuple to a complex-frequency input tuple
#	call as
#		tupC = jHWavesRtoC(tupR)
#
#	K. Brink 3/15/2026
 
	println(' ')
	println("Real frequncy was ", tupR[3], " rad/s")	
	wg = inputf(" Enter first guess at real and complex parts of frequency (2-elemnent array) (wg) (1/sec) ")

	ilong = NaN;

	tupC =(tupR[1:2]..., wg, tupR[4:6]..., ilong, tupR[8:end]...)
	
	return tupC

	end							# end of function






"""
	jHWavesCtoR(tupC)

	Convert the input tuple for complex frequency calculations
		to one for real frequency
	
K. Brink 3/15/2026
"""

function jHWavesCtoR(tupC)

#	Convert complex-frequency input tuple to a real-frequency input tuple
#	call as
#		tupR = jHWavesCtoR(tupC)
#
#	K. Brink 3/15/2026
 
	println(' ')
	wold = tupC[3] 
	println("Complex frequency was ", wold[1], "  ", wold[2], "im rad/s")	
	wg = inputf(" Enter first guess at the real frequency (wg) (1/sec) ")

	println(" ") 
	ilong = inputi(" Enter 1 for general frequency, wavenumber, 0 for long wave limit ")
		if abs(tupC[24]) > 0.001 && ilong == 0
			println(" ")
			println("Error! Can not use the long wave limit with a mean flow!")
			return NaN	
		end


	tupR =(tupC[1:2]..., wg, tupC[4:6]..., ilong, tupC[8:end]...)
	
	return tupR

	end							# end of function







"""
	jHWavesCSetup()

	Create the input tuple for use with jHSWavesCM

K.H. Brink 3/9/2026
"""
j

function jHWavesCSetup()

#   Create an input tuple for "jHWavesCM", given inputs
#   Call as
#       arr = jHWavesCSetup()
#   In response to queries, you can enter either numbers or array names
#
#   K.H. Brink  11/19/2025  based on Matlab 6/11/2018 version, 3/9/2026

	println(' ')
	println("      Coastal-Trapped waves with complex frequency")
	println(' ')
	println("This function will ask you a sequence of questions ")
	println( "     that are used to build the input tuple.     ")
	println(' ')
    

	nn = inputi("How many total gridpoints do you want in the cross shelf direction? (nn) ")
	if nn <5
    		println(' ')
    		println("Error! Must have nn > 4")
		return NaN
	end

	mm = inputi("How many total gridpoints do you want in the vertical? (mm) ")
	if mm < 5
    		println(' ')
    		println("Error! Must have mm >4")
		return NaN
	end
	arr = (nn, mm)

	wg = inputf("Enter first guess at real and complex parts of frequency (2-elemnent array) (wg) (1/sec) ")
		if length(wg) != 2
			println("Error: frequency should be an array of length 2!")
			return NaN
		end
	del = inputi("Enter 0 for a rigid lid, 1 for a free surface (del) ")
	icbc = inputi("Enter 0 for a closed x= 0 boundary, 1 for open (icbc) ")
	iobc = inputi("Enter 0 for a closed x =xmax boundary, 1 for an open (iobc) ")
	#ilw = inputi("Enter 1 for general frequency, wavenumber or 0 for long wave limit ")
	ilw = NaN
	f = inputf("Enter the Coriolis parameter (f) (rad/sec) ")
	xmax = inputf("Enter the domain width (xmax) (km) ")
	eps = inputf("Enter the nominal fractional accuracy for the frequency solution (eps) ")
	npts = inputi("Enter the number of frequencies to be computed (npts) ")

	#if ilw < 0.5
    	#       It only makes sense to compute one point for the long wave limit
    	#	npts = 1
	#end
	
	arr = (arr..., wg, del, icbc, iobc, ilw, f, xmax, eps, npts)					# 11 items



	rlz = inputf("Enter the first alongshore wavenumber to use (rlz) (rad/cm) ")
	if npts < 1.1
    		drl = rlz
	else
    		drl = inputf("Enter the wavenumber increment to use after rlz (drl) (rad/cm) ")
	end
	arr = (arr..., rlz,drl)										# 13 items

#   Read in depth
	ndep = inputi("How many distance, depth pairs will you provide (ndep >=1) ")
	
	xdep = inputf("Array of offshore distances for depth values (xdep in km) (dimension ndep) ")
	if ndep != length(xdep)
		println(' ')
		println("Error! array length must equal ndep!")
		return NaN
	end
	depr = inputf("Array of depths corresponding to xdep (depr in m) ")
	if ndep != length(depr)
		println(' ')
		println("Error! array length must equal ndep!")
		return NaN
	end
	if xdep[1] > 0.0001
		ndep = ndep + 1
		xdep = [0.0 xdep]
		depr = [depr[1] depr]
	end
	if xdep[end] < xmax
		ndep = ndep + 1
		xdep = [xdep xmax]
		depr = [depr depr[end]]
	end
	if ndep ==1
		ndep = 2
		xdep = [xdep xdep]
		depr = [depr depr]
	end
	iqq = findall(depr -> depr < 0, depr)
	if length(iqq) > 0.5
		println(' ')
		println("Error! All depth must be > 0!")
		return NaN
	end
	arr = (arr..., ndep, xdep, depr)							# 16 items
		


#       Read in bottom friction
	nr = inputi("Number of distance, bottom friction pairs to read (nr) ")
	if nr == 0
		nr = 2
		xr = [0 xmax]
		rr = [0. 0.]
	else
    		xr = inputf("Offshore distances for bottom friction values (xr in km) ")
    		rr = inputf("Array of bottom friction values corresponding to xr (rr in cm/sec) ")
	        if length(xr) != nr
        		println(' ')
        		println("Error! Array size must match nr!")
        		return NaN
    		end
    		if length(rr) != nr
        		println(' ')
        		println("Error! Array size must match nr!")
        		return NaN
	    	end
		if xr[1] > 0.0001
			nr = nr + 1
			xr = [0.0 xr]
			rr = [rr[1] rr]
		end
		if xr[end] < xmax
			nr = nr + 1
			xr = [xr xmax]
			rr = [rr rr[end]]
		end
									# 19 items
	end
    	arr = (arr..., nr, xr, rr)
  

#       Read in base-state stratification
	zr = inputf("Depth increment for input of Nsquared values? (zr in m) ")
	alph = inputf("Exponential tail length for Nsquared extrapolation (alph in km) ")
	nsqr = inputf("Nsquared values starting at the surface (nsqr in rad^2/sec^2) (nnsq values) ")
	nnsq = length(nsqr)

	arr = (arr..., nnsq, zr, alph, nsqr)							# 23 items
	
#       Read in mean flow

	vzero = inputf("Input peak value of mean alongshore flow (vzero: cm/sec) ")


	if abs(vzero) > 0.001
    		xzero = inputf(" Input distance offshore to peak mean flow (km) ")
    		zzero = inputf(" Input depth of peak mean flow (m) ")
    		zscaled = inputf(" Downward exponential scale of mean flow? (m) ")
    		zscaleup = inputf(" Upward exponential scale of mean flow? (m) ")
    		xscaleoff = inputf(" Offshore exponential scale of mean flow? (km) ")
    		xscaleon = inputf(" Onshore exponential scale of mean flow? (km) ")
    		kk = inputi(" Enter 1 for undisturbed Nsquared offshore, 0 for onshore ")
		nv = 1
	else
		xzero = 20.
		zzero = 0.
		zscaled = 10000.
		zscaleup = 10000.
		xscaleoff = 1000.
		xscaleon = 1000.
		kk = 1
		nv = 0
	end

    	arr = (arr..., vzero, xzero, zzero, zscaled, zscaleup, xscaleoff, xscaleon, kk)		# 31 items

#	Miscellaneous
	ipause = inputi("Enter 0 to skip pauses to see graphics, 1 to see graphics during execution ")
	arr = (arr..., ipause)									# 32 tiems total



#	This completes the reading in of information
#		Now for some graphics

	
	fig = Figure()
	display(GLMakie.Screen(),fig)

	colsize!(fig.layout, 1, Aspect(1, 1.8))
				
	ax1 = Axis(fig[1,1], title = "rplot")				# top panel
	ax2 = Axis(fig[2,1], title = "vplot")				# lower

	lines!(ax1,vec(xr),vec(rr))
	zc = maximum(rr)
	xlims!(ax1,0,xmax)
	if zc < 1e-8
		ylims!(ax1,0,0.1)
		text!(ax1,xmax/3,0.03, text = "r = 0")
	else
		ymaxx = 1.1*maximum(rr)
		ylims!(ax1,0, ymaxx)
	end
	hidedecorations!(ax1, grid = true, ticklabels = false, label = false)
	ax1.ylabel = "r [cm/s]"
	ax1.title = "Friction Coefficient"

	zmax = maximum(depr)
	yupper = -depr
	nzz = length(yupper)
	ylower =  -zmax*ones(nzz)
	band!(ax2,vec(xdep),vec(ylower),vec(yupper),color = (:blue,0.7))
	lines!(ax2,vec(xdep),vec(-depr), color = (:blue,0.7))
	xlims!(ax2,0,xmax)
	ylims!(ax2,-zmax,0)
	hidedecorations!(ax2, grid = true, ticklabels = false, label = false)
 	ax2.ylabel = "z [m]"
	ax2.xlabel = "x [km]"
	vstr = @sprintf("Depth and Mean v (vzero = %.1f cm/sec)", vzero )
	ax2.title = vstr

	if nv == 0
		
		text!(ax2,2*xmax/3,-0.5*zmax, text = "Mean v = 0")
	else
		nnn = 41
	    	mmm = 31
	        x = 0:(xmax/(nnn-1)):xmax
    	        z = -zmax:(zmax/(mmm-1)):0
    		zzero = -zzero

    		xfactoff = xscaleoff^2
    		xfacton = xscaleon^2
    		zfactu = zscaleup^2
    		zfactd = zscaled^2
    
    		vtemp = NaN*ones(mmm,nnn)
    
    		ii = findall(x-> x < xzero,x)
    		iic = findall(x -> x >=  xzero,x)
    		jj = findall(z -> z < zzero,z)
    		jjc = findall(z -> z >= zzero,z)
 
    		if isempty(size(ii)) != 1
        	   for n = ii[1]:maximum(ii)
			
            		if isempty(size(jj)) != 1
				sjj = length(jj)
				zzz = zzero*ones(sjj)
               			 vtemp[jj,n] = (exp(-((x[n]-xzero).*(x[n]-xzero)/xfacton)))*(exp.(-(z[jj]-zzz).*(z[jj]-zzz)/zfactd))
            		end
           	 	if isempty(size(jjc)) != 1
				sjj = length(jjc)
				zzz = zzero*ones(sjj)
                		vtemp[jjc,n] = (exp(-((x[n]-xzero).*(x[n]-xzero)/xfacton)))*(exp.(-(z[jjc]-zzz).*(z[jjc]-zzz)/zfactu))
            		end
        	   end
    		end
    
    		if isempty(size(iic)) != 1
        	   for n = minimum(iic):nnn
            		if isempty(size(jj)) != 1
				sjj = length(jj)
				zzz = zzero*ones(sjj)
                 		vtemp[jj,n] = (exp(-((x[n]-xzero).*(x[n]-xzero)/xfactoff)))*(exp.(-(z[jj]-zzz).*(z[jj]-zzz)/zfactd))
            		end
            		if isempty(size(jjc)) != 1
				sjj = length(jjc)
				zzz = zzero*ones(sjj)
                 		vtemp[jjc,n] = (exp(-((x[n]-xzero).*(x[n]-xzero)/xfactoff)))*(exp.(-(z[jjc]-zzz).*(z[jjc]-zzz)/zfactu))
            		end
        	   end
   	 	end
		vtemp = vzero*vtemp
    
    
    		contour!(ax2,x, z, vtemp', levels = 5, color = :black, labels = false)
    		band!(ax2,vec(xdep),vec(ylower),vec(yupper),color = (:blue,0.7))
	end
	display(fig)

	vtemp = nothing
	return arr
end									# end of function



"""
	jHWavesCFinch(arr)

	Make changes to input tuple for complex frequency
	Returns a revised tuple

K.H. Brink 11/21/2025
"""


function  jHWavesCFinch(arr)

#   Function to modify an existing jHWavesC input file conveniently
#
#   K.H. Brink,  3/9/2026, 4/2/2026

	println(' ')
	println("     Function to create a modified input tuple")
	println(" ")

	xx = """First you need to select what you want to change.
                                                      
          Options are:                                
            g:  Grid size                     
            w:  Initial Frequency guess                                 
            b:  Boundary conditions
            f:  Coriolis parameter                     
            x:  Domain size                        
            d:  Dispersion curve definition
            e:  Nominal accuracy               
            h:  Depth profile                          
            r:  Bottom friction
            n:  Stratification
            v:  Mean flow field
            p:  Pauses during execution
                                                      
     Any arrays are row arrays, not column arrays   

     Please type in the appropriate letter and hit return  """  
                                                      
 	println(xx)

	stinn = readline()
	stin = stinn[1]
	println(' ')   

	if stin == 'g'
    		println("Old [nn   mm] = ", arr[1], ",  ", arr[2])
    		iii = inputi(" Enter new [nn  mm] ")
		nn = iii[1]
		mm = iii[2]
    		narr = ( nn, mm, arr[3:end]...)

	elseif stin == 'w'
    		println("Old frequency guess (real and complex parts) (sec^-1) ",  arr[3])
    		www = inputf(" Enter new real and complex frequency parts (2-element array) (sec^-1) ")
		if length(www) != 2
			println(' ')
			println("Error: frequency should be an array of length 2!")
			return
		end
		narr = (arr[1:2]..., www, arr[4:end]...)
    

	elseif stin == 'b'
    		println("Old boundary conditions")
    		if arr[4] == 0
        		println("   Rigid lid")
    		else
        		println("   Free surface")
    		end
    		if arr[5] == 0
        		println("   Solid coastal wall")
    		else
        		println("   Open coastal boundary")
    		end
    		if arr[6] == 0
        		println("   Closed boundary at x = xmax")
    		else
        		println("   Open boundary at x = xmax")
    		end
    		ar1 = inputi(" Enter 0 for rigid lid, 1 for free surface ")
   		ar2 = inputi(" Enter 0 for no flow through coastal wall, 1 for open boundary ")
    		ar3 = inputi(" Enter 0 for closed offshore boundary, 1 for open ")
    		narr = (arr[1:3]..., ar1, ar2, ar3, arr[7:end]...)
    
    		if ar2 > 0.5
        		println("Note: bottom must be flat at x = 0")
    		end
		if ar3 > 0.5
			println("Note: bottom must be flat at x = xmax")
		end

	elseif stin == 'f'
    		println("Old Coriolis parameter = ", arr[8])
    		www = inputf(" Enter new Coriolis parameter  (sec^-1) ")
		narr = (arr[1:7]..., www, arr[9:end]...)

	elseif stin == 'x'
	    	println("Old offshore domain width (km) = ", arr[9])
    		www = inputf(" Enter new domain size (km) ")
		narr = (arr[1:8]..., www, arr[10:end]...)
		

	elseif stin == 'e'
    		println("Old nominal fractional accuracy for frequency= ", arr[10])
    		www = inputf(" Enter new nominal fractional accuracy ")
		narr = (arr[1:9]..., www, arr[11:end]...)

	elseif stin == 'd'
    		println("Previous number of dispersion curve points = ", arr[11])
    		iii = inputi(" Enter the new number of points on the dispersion curve ")
    		println("Old first wavenumber, increment along curve (cm^-1) = ", arr[12:13])
    		www = inputf(" Enter new first alongshore wavelength (cm-1) ")
    		if iii < 1.1
        		ar1 = arr[13]
    		else
        		ar1 = inputf(" Enter new wavenumber increment (cm-1) ")
    		end
		if arr[7] < 0.5								
			iii = 1
			println("Only one wavelength is used in the coastal long wave approximation")
		end

		narr = (arr[1:10]..., iii, www, ar1, arr[14:end]...)

	elseif stin == 'h'
    		println("Old number of input points = ", arr[14])
    		nd = inputi(" Enter number of (x, h) pairs  (> 0) ")
    		if nd < 1
			println(' ')
        		println("Error!! nd must be >= 1")
			return
    		end
    		println("Old array of x locations for depth (km) = ", arr[15])	
    		ar1 = inputf(" Enter array of x locations where depth will be given (km) ")

    		println("Old array of depths (m) = ", arr[16])
    		ar2 = inputf(" Enter array of depths at locations given (m) ")
		
		iqq = findall(ar2 -> ar2 < 0, ar2)
		if length(iqq) > 0.5
			println(' ')
			println("Error! All depth must be > 0!")
			return
		end



		if ar1[1] > 0.
			nd = nd + 1
			ar1 = [0.0 ar1]
			ar2 = [ar2[1] ar2]
		end
		if ar1[end] < arr[9]
			nd = nd + 1
			ar1 = [ar1 arr[9]]
			ar2 = [ar2 ar2[end]]
		end
		narr = (arr[1:13]..., nd, ar1, ar2, arr[17:end]...)


	elseif  stin == 'r'

    		println("Old number of r values = ", arr[17])
    		nr = inputi(" Enter number of (x, r) pairs (>=0) ")
  		if nr !=0
        	        println("Old r locations = ", arr[18])
        		ar1 = inputf(" Enter array of x locations where resistance coefficient is given (km) ")
            		println("Old r values = ", arr[19])
			ar2 = inputf(" Enter array of resistance coefficient (cm/sec) ")
			if length(ar1)!= nr
				println(' ')
				println("Error! Array length not equal to nr!")
				return
			end
			if length(ar2)!= nr
				println(' ')
				println("Error! Array length not equal to nr")
				return
			end
        	else
			nr = 2
			ar1 = [0 arr[9]]
			ar2 = [0 0]
		end
		narr = (arr[1:16]..., nr, ar1, ar2, arr[20:end]...)		
		
       
	elseif stin =='n'
		println("Old # of N^2 points = ", arr[20]) 
      		nn = inputi(" Enter new number of Nsquared values to read (>=1) ")
		println("Old depth increment for reading N^2 (m) = ", arr[21])
		ar1 = inputf(" Enter new depth increment for reading N^2 (m) ")
		println("Old exponential tail for need N^2 profile (km) = ", arr[22])
		ar2 = inputf(" Enter new exponential tail for N^2 (km) ")
		println("Old array of N^2 values (sec^-2) = ", arr[23])
		ar3 = inputf(" Enter new array of N^2 values (sec^-2) ")
		nn = length(ar3)
		if length(ar3) != nn
			println(' ')
			println("Error! Array length not equal to nn!")
			return
		end
		narr = (arr[1:19]..., nn, ar1, ar2, ar3, arr[24:end]...)
    
	elseif stin == 'v'

		println("Old alongshore flow amplitude (cm/sec) = ", arr[24])	
 		va = inputf(" Enter mean alongshore flow amplitude (cm/sec) ")
		if abs(va) > 0.0001
			println("Old offshore distance of speed maximum (km) ", arr[25])
		        ar1 = inputf(" Enter new offshore distance of velocity extreme (km) " )   
			println("Old vertical location of velocity extreme (m) = ", arr[26])
			ar2 = inputf(" Enter new vertical location of velocity extreme (m) ")
			println("Old downward and upward scales for mean velocity (m) ", arr[27], ",   ", arr[28])
			ar3 = inputf(" Enter new downward scale for mean flow (m) ")
			ar4 = inputf(" Enter new upward scale for mean flow ")
			println("Old offshore and onshore scales for mean flow (km) ", arr[29], ",  ", arr[30])
			ar5 = inputf(" Enter new offshore scale for mean velocity (km) ")
			ar6 = inputf(" Enter new onshore scale for mean velocity (km) ")
			if arr[31] < 0.5
				println("Old: Stratification is undisturbed at the coast")
			else
				println("Old: Stratification is undisturbed far offshore")
			end
			istrat = inputi(" Enter 0 for stratification undisturbed at the coast, 1 for undisturbed offshore ")
		else
			ar1 = 30.
			ar2 = 0.
			ar3 = 1000.
			ar4 = 1000.
			ar5 = 1000.
			ar6 = 1000.
			istrat = 1.
		end
		narr = (arr[1:23]..., va, ar1, ar2, ar3, ar4, ar5, ar6, istrat, arr[32]...)
    
	elseif stin =='p'
    	
    		if arr[end] > 0.5
		        println("Was set to pause for graphics")
		else
        		println("Was set not to pause for graphics")
    		end
    		ipause = inputi(" Enter 1 to pause for graphics or 0 to skip pauses ")

    		narr = (arr[1:31]..., ipause)
    
	else	
		println(' ')
    		println("Value entered does not correspond to any of the above categories")
		return
	end
	return narr

end									# end of function








"""
    inputi(quer)
    
    Replace the Matlab "input" function for integer numbers/arrays
  
  	quar = a string that asks for numbers/arrays

  It asks for numbers. You can either:
	type in a single number
	type in an array, e.g.:  [1 2 3]
     or (if you have the numbers in an array or vector x)
	type in string(x)
	
  
  Returns an integer number or an array of integers
  
  K. Brink  10/21/2025
"""

function inputi(quer)
	local respp, iar, icount, contt


    	print(quer)
    										#respp is a string
    	respp = readline()

	if length(respp) >= 6
		if respp[1:6] == "string"					# This allows you to enter an array by name
			ccc = Meta.parse(respp)
			rrr = eval(ccc)
			respp = rrr
		end
	end
    
    	while respp[1] == ' '            					# Remove any initial blanks
        	respp = respp[2:end]
    	end
    	while respp[end] == ' '	 						# remove any blanks at the end
		respp = respp[1:(end-1)]
	end
	
	
	if occursin(",",respp)							# remove commas and replace with blanks
		
		respp = replace(respp,',' => ' ')
	end
	
    
    	if respp[1] != '['							# single number, no brackets
        	iar = parse(Int,respp)
        	return iar
    	else									# input is a vector, not a scalar
		if respp[end] != ']'
			println("Error: missing bracket!")
			return NaN
		else
			respp = respp[2:(end-1)]    		 		#remove brackets
			while respp[end] ==' '
				respp = respp[1:(end-1)]
			end

			icount = 1;
			contt = 1;
			while contt > 0.5
				if occursin(' ',respp)
					ii = findfirst(' ',respp)
					tem = respp[1:(ii-1)];
					respp = respp[(ii+1):end]		# Shorten string
					while respp[1] == ' '
						respp = respp[2:end]
				 	end
				 	iii = parse(Int64,tem)
				 	if icount == 1
						iar = iii
				 	else
						iar = [iar iii]
				 	end
					icount = icount + 1
				else						# You are at the last number
					contt = 0				# signal that you are done
					iii = parse(Int64,respp)
					if icount > 1.5
						iar = [iar iii]
					else
						iar = iii
					end
				end
			end					
		
		return iar
		end						
	end							
	
end										# function end




"""
    inputf(quer)
    
  Replace the Matlab "input" function for floating point numbers/arrays
  
  	quer = a string that asks for numbers/arrays, e.g., 'Array of x (km)? "

  It asks for numbers. You can either:
	type in a single number
	type in an array, e.g.:  [1.1 2.6 3.2]
     or (if you have the numbers in an array or vector x)
	type in string(x)
  
  Returns a number or an array of floating point numbers
  
  K. Brink  10/21/2025
"""

function inputf(quer)
	local respp, flar, icount, contt
    	print(quer)
    
   	respp = readline()						# this is a string

	if length(respp) >= 6
		if respp[1:6] == "string"				# this allows you to enter an array by name
			ccc = Meta.parse(respp)
			rrr = eval(ccc)
			respp = rrr
		end
	end
    
    	while respp[1] == ' '            				# Remove any initial blanks
        	respp = respp[2:end]
    	end
	while respp[end] == ' '						# remove any blanks at the end
		respp = respp[1:(end-1)]
	end
	
	if occursin(",",respp)						# remove commas and replace with blanks
		
		respp = replace(respp,',' => ' ')
	end
	
		
	
    
    	if respp[1] != '['						# single number, no brackets
        	flar = parse(Float64,respp)
        	return flar
    	else								# input is a vector, not a scalar
		if respp[end] != ']'
			println("Error: missing bracket!")
			return NaN
		else
			respp = respp[2:(end-1)]   	 	 	#remove brackets
			while respp[end] ==' '
				respp = respp[1:(end-1)]
			end

			icount = 1
			contt = 1
			while contt > 0.5
				if occursin(' ',respp)
					ii = findfirst(' ',respp)
					tem = respp[1:(ii-1)]
					respp = respp[(ii+1):end]		# Shorten string
					while respp[1] == ' '
						respp = respp[2:end]
				    	end
					fii = parse(Float64,tem)
					if icount == 1
						flar = fii
					else
						flar = [flar fii]
					end
					icount = icount + 1
				else						# You are at the last number
					contt = 0				# signal that you are done
					fii = parse(Float64,respp)
					if icount > 1.5
						flar = [flar fii]
					else
						flar = fii
					end
				end
			end
			return flar
		end
	end
end										# function end



