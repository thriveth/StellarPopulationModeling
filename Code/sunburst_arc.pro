function gconv, x, sigma, edge_wrap=edge_wrap, fwhm = fwhm
  if (n_elements(fwhm) gt 0) then $
    sigma = fwhm/2.3548
  if sigma eq 0 then return, x
  binfactor=1
  ksize=round(4.0*sigma+1.0)*2
  xx = findgen(ksize)-ksize/2
  kernel=exp(-xx^2/(2*sigma^2))
  kernel=kernel/total(kernel)
  sm = convol(x, kernel, edge_wrap=edge_wrap)
  return, sm
end



function reddy_unred, wave, flux, ebv, funred, R_V = R_V
  x  = wave*0.0001    
  klam = 2.191 +.974/x
  funred = flux*10.0^(0.4*klam*ebv)
end



function s99_mcombine, x, a, mlib=mlib
  y = double(mlib) # a[1:*]
  ; Calzetti reddening
  ;calz_unred, x, y, -1.*a[0], ynew
  ; Reddy reddening
  reddy_unred, x, y, -1.*a[0], ynew
  ;Using an LMC from Gordon et al. 2003
  ;fm_unred, x, y, -a[0], ynew, c1 = -4.959, c2 = 2.264, c3 = 0.389, c4=.461, x= 4.6, gamma = 1
  return, ynew

end


function s99_continuum, model, restwl, flux, err, vdisp, $
  yfit = yfit, noplot=noplot, struct_only = struct_only

  nage = n_elements(model.age)

  coefs = {ebv: 0.0, ebv_err: 0.0, $
    light_frac: fltarr(nage), light_frac_err: fltarr(nage), $
    model_age: model.age, vdisp: 0.0, vdisp_err: 0.}

  nmodels = n_elements(model.age)

  parinfo = replicate({value:0.D, fixed:0, limited:[0,0], tied:'', $
    limits:[0.D,0.D]}, nmodels+1)

  parinfo[0].limited = [1,1]
  parinfo[0].fixed =0
  parinfo[0].limits = [0.0,5.0] 
  parinfo[0].value = 0.3
  parinfo[1:*].limited = [1,0] 
  parinfo[1:*].limits = [0.0,0.0]
  parinfo[1:*].value = 0.1

  npix = n_elements(restwl)
  quality = fltarr(npix) + 1

  outside_model = where(restwl le min(model.wave) or restwl ge max(model.wave))
  if outside_model[0] ne -1 then quality[outside_model] = 0

  bad = where(finite(flux) ne 1 or finite(err) ne 1 or err lt 1e-6 or err eq 0)
  if bad[0] ne -1 then quality[bad] = 0

  ok = where(quality eq 1)
  s99_pix = 146.107 
  s99_vdisp = 92.0 

  if vdisp lt s99_vdisp then vdisp_add = 0 $
  else vdisp_add = sqrt(vdisp^2 - s99_vdisp^2)
  sigma_pix = vdisp_add / s99_pix

  custom_lib = dblarr(npix, nmodels)
  for ii = 0, nmodels - 1 do begin
    cflux = gconv(model.flux[*,ii], sigma_pix) 
    custom_lib[*,ii] = interpol(cflux, model.wave, restwl)
    custom_lib[*, ii] = custom_lib[*, ii]/median(custom_lib[where(restwl gt 1267 and restwl lt 1276), ii])
 endfor
  if outside_model[0] ne -1 then custom_lib[outside_model, *] = 0.0

  fitcoefs = mpfitfun('s99_mcombine', restwl[ok], flux[ok], err[ok], $
    parinfo = parinfo, functargs = {mlib: custom_lib[ok,*]}, $
    perror=perror, niter=niter, status=status, /quiet, $
    maxiter = 1000, errmsg = errmsg)
  yfit = s99_mcombine(restwl, fitcoefs, mlib=custom_lib)
  return, coefs

end

pro sunburst_arc


cspeed = 2.9979E5

cd, current = cur_dir

!PATH = !PATH +':'+Expand_path('+'+cur_dir) + Expand_path('+' + cur_dir + '/MPFIT/)
gal = 'Sunburst_Arc'

z = 2.37094
  
readcol, cur_dir + 'Sunburst_Arc.txt', owave1, iflux2, ierr2

wave = owave1 / (1.+z)

iflux1 = iflux2 *(2.99E18)/owave1^2
ierr1 = ierr2 * (2.99E18)/owave1^2

nflux1 = iflux1/median(iflux1[where(wave gt 1267 and wave lt 1276)])
nerr1 = ierr1/median(iflux1[where(wave gt 1267 and wave lt 1276)])
super_mods= mrdfits(cur_dir+'super_mods.fits', 1)

s1 = {wave: wave, flux: nflux1, err: nerr1}
s = create_struct(s1, 'name', gal, 'z', z)
readcol, cur_dir+'Sunburst_Arc.linelist', $
  mwlines1, id, dumbywave, fval, dumby, color, zoff, type, format = 'f, a, f, f,f, a, f, a'
dv = 500.
mask = fltarr(n_elements(wave))+1.
mwlines = mwlines1[where(strmatch(type, 'ISM') eq 1 or strmatch(type, 'INTERVE') eq 1 or strmatch(type, 'EMISSION') eq 1)]

mwlines = [mwlines, 1375, 1377, 1380, 1407, 1465, 1467, 1535, 1545, 1688, 1747, 1867, 1876, 1962]

nmwlines = n_elements(mwlines)

for maskn = 0, nmwlines-1 do begin & $
  tempmask = where(wave gt -dv/cspeed*mwlines[maskn]+mwlines[maskn] and wave lt dv/cspeed*mwlines[maskn]+mwlines[maskn])   & $
  if mwlines[maskn] eq 1550.77 then tempmask = where(wave gt -dv/cspeed*mwlines[maskn]+mwlines[maskn] and wave lt 50./cspeed*mwlines[maskn]+mwlines[maskn])   & $
mask[tempmask] = 0 & $
endfor

chisreg = where(mask eq 1  and s.wave gt 1240 and s.wave lt 1900 and finite(s.flux) eq 1 and s.err ne 0, comp = badc)
s.err[badc] = 0

vdisp = 151.01721

c_super_mods = s99_continuum(super_mods, s.wave, s.flux, s.err, vdisp, yfit = cfit_super)
exz = total(c_super_mods.light_frac *super_mods.z )/total(c_super_mods.light_Frac )
exage = total(c_super_mods.light_frac *super_mods.age )/total(c_super_mods.light_Frac )

end

