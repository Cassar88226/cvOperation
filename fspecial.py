import numpy as np
import math
# define fucntion fspecial(shape = disk)
def fspecial(type, radius):
    if type == 'disk':
        rad   = radius
        crad  = math.ceil(rad-0.5)
        x = np.arange(-crad, crad + 1)
        y = np.arange(-crad, crad + 1)
        xx,yy = np.meshgrid(x, y, sparse = True)
        


        maxxy = np.maximum(np.absolute(xx),np.absolute(yy))
        minxy = np.minimum(np.absolute(xx),np.absolute(yy))

        a = (rad**2 <  (maxxy+0.5)**2 + (minxy-0.5)**2) 
        b = (rad**2 >= (maxxy+0.5)**2 + (minxy-0.5)**2)    
        m1 = (a * 1 * (minxy-0.5) + b * 1 *np.sqrt((rad**2 - np.square(maxxy + 0.5)).astype(complex))).real
        
        a1 = (rad**2 >  (maxxy-0.5)**2 + (minxy+0.5)**2)
        b1 = (rad**2 <= (maxxy-0.5)**2 + (minxy+0.5)**2)    
        m2 = (a1 * 1 * (minxy+0.5) + b1 * 1 * np.sqrt((rad**2 - np.square(maxxy - 0.5)).astype(complex))).real

        
        a2 = (rad**2 <  (maxxy+0.5)**2 + (minxy+0.5)**2) 
        b2 = (rad**2 >= (maxxy-0.5)**2 + (minxy-0.5)**2)
    #     sgrid = (rad^2*(0.5*(asin(m2/rad) - asin(m1/rad)) + ...
    #             0.25*(sin(2*asin(m2/rad)) - sin(2*asin(m1/rad)))) - ...
    #             (maxxy-0.5).*(m2-m1) + (m1-minxy+0.5)) ...
    #             .*((((rad^2 < (maxxy+0.5).^2 + (minxy+0.5).^2) & ...
    #             (rad^2 > (maxxy-0.5).^2 + (minxy-0.5).^2)) | ...
    #             ((minxy==0)&(maxxy-0.5 < rad)&(maxxy+0.5>=rad))));
        f1 = (0.5*(np.arcsin(m2/rad) - np.arcsin(m1/rad)) + 0.25*(np.sin(2*np.arcsin(m2/rad)) - np.sin(2*np.arcsin(m1/rad))))
        first = (rad**2*f1 - (maxxy-0.5)*(m2-m1) + (m1-minxy+0.5))
        second = (((a2 & b2) | ((minxy==0)&(maxxy-0.5 < rad)&(maxxy+0.5>=rad))))

        sgrid = second * 1 * first
        
        a3 = (rad**2 >  (maxxy+0.5)**2 + (minxy+0.5)**2)

        sgrid = sgrid + a3 * 1    
        sgrid[crad,crad] = min(np.pi*rad**2,np.pi/2);
        
        
        if ((crad>0) and (rad > crad-0.5) and (rad**2 < (crad-0.5)**2+0.25)):
            
            m1  = sqrt(rad**2 - (crad - 0.5)**2)
            m1n = m1/rad
            sg0 = 2*(rad**2*(0.5*np.arcsin(m1n) + 0.25*np.sin(2*np.arcsin(m1n)))-m1*(crad-0.5))
            sgrid[2*crad,crad]     = sg0
            sgrid[crad,2*crad]     = sg0
            sgrid[crad,0]          = sg0
            sgrid[0,crad]          = sg0
            sgrid[2*crad-1,crad]   = sgrid(2*crad-1,crad) - sg0
            sgrid[crad,2*crad-1]   = sgrid(crad,2*crad-1) - sg0
            sgrid[crad,1]          = sgrid(crad,1)        - sg0
            sgrid[1,crad]          = sgrid(1,crad)        - sg0
        
        
        sgrid[crad,crad] = min(sgrid[crad,crad],1);
        
        h = sgrid/np.sum(sgrid);
    return h


    