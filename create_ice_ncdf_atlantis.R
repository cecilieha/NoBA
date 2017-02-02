#NB!!!!!!!!!
# after running this program you need to rename the timedimension, as R doesn't allow
# same names of dimensions and variables. Use nco and write the following
# ncrename -d t1,t filename.nc newfilename.nc, and you're ready to go.
#libraries
library(ncdf4)
library(pracma)
library(maptools)
if (!rgeosStatus()) gpclibPermit()
library(methods)
library(sp)


# load in the polygon (lon, lat) corner coordinate locations:
cb <- readShapePoly("/data/NoBaAtlantis/Common_files/Grid/MENUIIareasPolNewId_grass_tol0p01.shp")
box_id<-sort(as.numeric(paste(unique(cb$box_id))))                                                        ########
Nt_get = c(1)		# Time ranges to process from NetCDF file - ant tidssteg

year_all=c(2006:2025) #define over which period you'd like to calculate the forcing files (each year will be written in a separate file)

corners <- read.table("/data/NoBaAtlantis/Common_files/corners_neighbours_nordic.txt", header=T)
boxes = box_id+1 #nr of polygons

#at most 6 timesteps in each file, hence
# calculate on subset instead of whole area 
x_grab=c(380,770)
y_grab=c(570,950)
#NB! Model grid dependent
x1=1
x2=390
xt=x2-x1
  #
y1=1  
y2=380
yt=y2-y1 

#the info on landmask, csr etc should be part of your romsfile - check
load('/data/NoBaAtlantis/Common_files/missing_roms_vars.RData')

mask.file="/data/NoBaAtlantis/Common_files/AA_10km_grid.nc"
m.f=nc_open(mask.file)
#mask_rho = ncvar_get(m.f, "mask_rho") # mask on RHO-points; all 0 or 1, dim( 507 329 )
lon_rho  = ncvar_get(m.f, "lon_rho")
lat_rho  = ncvar_get(m.f, "lat_rho")
h        = ncvar_get(m.f, "h" ) # 
nc_close(m.f)

for (y in year_all){ #NB! valid if you have forcing files with year included in names, otherwise rewrite the loop 
  print(y)
  cnt=1
  f.path=paste("/data/NoBaAtlantis/Forcing/NorESM_2006_2070/NorESM_year_",y,'.nc',sep="")
  ex.cb = nc_open(f.path)   
   ocean_time = ncvar_get( ex.cb, "ocean_time" )  # 
   icea = ncvar_get( ex.cb, "aice" ) # fraction of cell covered by ice 
   iceh = ncvar_get( ex.cb, "hice" ) # ice thickness
   nc_close(ex.cb)

   # Get the horizontal, vertical, and time dimensions of the data
   Nt=length(ocean_time)

#  intialize arrays
   max_Ni=6
   t_icea=array(NA,dim=c(length(box_id),Nt*max_Ni))
   t_iceh=array(NA,dim=c(length(box_id),Nt*max_Ni))
   w=seq(0,1,by=(1/(ceil(365/Nt)))) #linear weights for converting from 5-day means to daily fields

   # if the time record chosen is not in the file, then exit:
   if (max(Nt_get) > Nt) {
      print('The time record exceeds the available number ....')	
   } else {
      print('Ready to Go!')
   }

   # Calculate the number of time intervals to average u, v, T and S over:
   Ti=Nt # Several time steps in each file
      

   # define arrays at internal rho-points to average velocities to and get the
   # variables at these point locations
   x_sub=array(1,dim=c(xt-2,yt-2))
   y_sub=array(1,dim=c(xt-2,yt-2))
   icea_sub=array(1,dim=c(xt-2,yt-2,Ti))
   iceh_sub=array(1,dim=c(xt-2,yt-2,Ti))

   icea_tmp=array(1,dim=c(xt-2,yt-2,Nt))  #need Nt as for some reason several days are stored in same file.
   iceh_tmp=array(1,dim=c(xt-2,yt-2,Nt))

   x_sub=lon_rho[(x_grab[1]+2):(x_grab[2]-2),(y_grab[1]+2):(y_grab[2]-2)]
   y_sub=lat_rho[(x_grab[1]+2):(x_grab[2]-2),(y_grab[1]+2):(y_grab[2]-2)]

   tempicea= icea 
   tempiceh= iceh 

   rm(icea,iceh)

   #Exclude land velocity values and set the .. to zero if necessary:

   II= which(abs(tempicea)>1000000) # NB!!! this can be model dependent! - check your fill values 
   if (length(II)!=0) tempicea[II]=0

   II= which(abs(tempiceh)>100000) #NB!!!!! this can be model dependent - check fill values!
   if (length(II)!=0) tempiceh[II] =0

   # Average temperature and salinities to internal rho-point locations:

   # For checking the values
   # filled.contour(temp_tmp[,,30,1]) #nb! roms is "upside-down", hence layer 30 is surface.

   icea_tmp=tempicea[(x1+1):(x2-2),(y1+1):(y2-2),]
   iceh_tmp=tempiceh[(x1+1):(x2-2),(y1+1):(y2-2),]

#remove NAs
   icea_tmp[is.na(icea_tmp)]=1e37
   iceh_tmp[is.na(iceh_tmp)]=1e37

###################


# Extract the finite velocity values and their (lon, lat) locations

   if (Ti == 1 ) TiLoop = 1
   if (Ti > 1 ) TiLoop = Ti 

   len = length(box_id)
   final = (array(0, dim=c(length(box_id),4,TiLoop)))


   for (tm in 1:TiLoop) { # loop by time
      qq = iceh_tmp[,,tm] # get the depth and time levels looping through now
      if (length(which(is.finite(qq)))>0) {

         II=is.finite(qq) # true false

         x_fin=x_sub[II] #  lon's
         y_fin=y_sub[II] # lats
         qq=drop(icea_tmp[,,tm]);icea_fin=qq[II]
         qq=drop(iceh_tmp[,,tm]);iceh_fin=qq[II]

         coordssat=SpatialPoints(coord=cbind(x_fin,y_fin))

         a=0
         for (i in 1:length(box_id)) { # 
	# this gives me the index of the cb section, not the JJ!
            cb_get<- which(cb$box_id == box_id[i]) # 
            a = a+1
            if (length(cb_get) > 0) {
	       JJ_get<-as.numeric(unlist(over(cb[cb_get,],coordssat,returnList=T))) # 

	       x_mean=mean(x_fin[JJ_get], na.rm=TRUE)
	       y_mean=mean(y_fin[JJ_get], na.rm=TRUE)
	       iceh_mean=mean(iceh_fin[JJ_get], na.rm=TRUE)
	       icea_mean=mean(icea_fin[JJ_get], na.rm=TRUE)

	       final[a,1,tm]=x_mean
	       final[a,2,tm]=y_mean
	       final[a,3,tm]=iceh_mean
	       final[a,4,tm]=icea_mean
	}

      }
   }
   } #end of time loop


   f_icea=drop(final[,4,])
   f_iceh=drop(final[,3,])

#####################   
#program to interpolate from 3-day means to daily values

   t_icea[,cnt:(cnt+Ti-1)]=f_icea
   t_iceh[,cnt:(cnt+Ti-1)]=f_iceh
   cnt=cnt+Ti
#} # loop by day   

#remove access entries
t_icea=t_icea[,1:(cnt-1)]
t_iceh=t_iceh[,1:(cnt-1)]

#
daymean=ceil(365/Nt)

#interpolate between - to get daily values. Each field is a mean over the number of days specified in daymean
icea_fin=array(NA,dim=c(length(box_id),(cnt-1)*daymean))
iceh_fin=array(NA,dim=c(length(box_id),(cnt-1)*daymean))


jcnt=1
for (i in 1:(cnt-2)){
for(j in 1:(length(w)-1)){
  icea_fin[,jcnt+j-1]=t_icea[,i]*w[(length(w)-j+1)]+t_icea[,i+1]*w[j]  
  iceh_fin[,jcnt+j-1]=t_iceh[,i]*w[(length(w)-j+1)]+t_iceh[,i+1]*w[j]  
}
  jcnt=jcnt+daymean
}

if(((jcnt-daymean))<365 & cnt>360){ #length other than 365 messes things up in long runs
   icea_fin[,(jcnt:365)]=icea_fin[,jcnt-cei]
   iceh_fin[,(jcnt:365)]=iceh_fin[,jcnt-daymean]
}

#remove fill-values
icea_fin[which(icea_fin>1)]=0.
iceh_fin[which(iceh_fin>1000)]=0.

#temporary arrays
icea_nc_temp=array(NA,dim=c((length(box_id)),dim(icea_fin)[2])) # max neighbours set to 20 for now
iceh_nc_temp=array(NA,dim=c((length(box_id)),dim(iceh_fin)[2])) # max neighbours set to 20 for now

#icea_nc_temp[,]=aperm(icea_fin,c(2,1))
#iceh_nc_temp[,]=aperm(iceh_fin,c(2,1))

icea_nc_temp=icea_fin
iceh_nc_temp=iceh_fin

#create time - easier than to store the info all the way. 
t_start=ocean_time[1]
dt=86400
t_stop=ocean_time[1]+(dim(icea_nc_temp)[2]-1)*86400
t_tot=seq(t_start,t_stop,dt)


#filename
f.ice=paste("NoBa_NorESM_ice_",y,".nc",sep="")

#define dimensions
dimb=ncdim_def("b","",1:length(box_id),create_dimvar=FALSE)
dimt=ncdim_def("t1","",1:length(t_tot),unlim=TRUE)#,create_dimvar=FALSE)

#create variables
#NB!!!!!! Unlimited rec needs to be on the right - otherwise the program complains!
origMissval=0.
var.t=ncvar_def("t","seconds since 2008-01-01 00:00:00 +10",dimt,0,prec="double")
var.icea=ncvar_def("Ice_Class1","frac",list(dimb,dimt),0,prec="double",origMissval)
var.iceh=ncvar_def("total_depth","m",list(dimb,dimt),0,prec="double",origMissval)

#create file
nc_ice=nc_create(f.ice,list(var.t,var.icea,var.iceh))

#assign global attributes to temp file
ncatt_put(nc_ice,0,"title","Temperature file, NoBa")
ncatt_put(nc_ice,0,"geometry","Nordic.bgm")
ncatt_put(nc_ice,0,"parameters","")

#assign attributes to variables
ncatt_put(nc_ice,var.t,"dt",86400,prec="double")

#assign variables to file
ncvar_put(nc_ice,var.icea,icea_nc_temp)
ncvar_put(nc_ice,var.iceh,iceh_nc_temp)

nc_close(nc_ice)
}
