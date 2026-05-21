import numpy

def wellen_simulieren(u, v, n_t, d_t, d_x, d_y):
    if numpy.max(v)*d_t/d_x > 1: #Sicherstellen, dass die Simulation möglich ist
        print(numpy.max(v)*d_t/d_x)
        print("v*d_t/d_x darf nicht größer als 1 sein. Sonst ist v schneller, als die Elemente nacheinander berechnet werden können.")
        
    for i in range(n_t-2):
        print(str(i+2)+" von "+str(n_t), end='\r')
        u[:,:, i+2] = future_wavefield_vek(v, u[:, :, i], u[:, :, i+1], d_t, d_x, d_y)
    return u

def future_wavefield(v, u_p, current_wavefield, d_t, d_x, d_y): #u_p = u(k,l,i-1), u_c = u(k,l,i)
    future_wavefield = numpy.zeros((n_x, n_y))
    u_c = numpy.zeros((n_x+2, n_y+2))
    u_c[1:-1,1:-1] = current_wavefield #zum differenzieren erweitern, damit das ergebnis die selben dimensionen hat
    for k in range(n_x):
        for l in range(n_y):
            x = ((numpy.square(v[k,l])*numpy.square(d_t))/numpy.square(d_x))*(u_c[k+1,l] - 2*u_c[k,l] + u_c[k-1,l])
            y = ((numpy.square(v[k,l])*numpy.square(d_t))/numpy.square(d_y))*(u_c[k,l+1] - 2*u_c[k,l] + u_c[k,l-1])
            future_wavefield[k,l] = 2*u_c[k,l] - u_p[k,l] + x + y            
    return future_wavefield

def future_wavefield_vek(v, u_p, c_w, d_t, d_x, d_y):   
    u_c = numpy.pad(c_w, ((1, 1), (1, 1)), mode='constant')  # Zero-padding

    x_diff = ((v * d_t / d_x) ** 2) * (u_c[2:, 1:-1] - 2 * u_c[1:-1, 1:-1] + u_c[:-2, 1:-1])
    y_diff = ((v * d_t / d_y) ** 2) * (u_c[1:-1, 2:] - 2 * u_c[1:-1, 1:-1] + u_c[1:-1, :-2])

    future_wavefield = 2 * u_c[1:-1, 1:-1] - u_p + x_diff + y_diff

    future_wavefield[:,-1] = c_w[:,-1] - (v[:,-2]*d_t/d_x)*(c_w[:,-1]-c_w[:,-2])
    #future_wavefield[0,:] = c_w[0,:] - (v[1,:]*d_t/d_x)*(c_w[0,:]-c_w[1,:])
    future_wavefield[:,0] = c_w[:,0] - (v[:,1]*d_t/d_x)*(c_w[:,0]-c_w[:,1])
    future_wavefield[-1,:] = c_w[-1,:] - (v[-2,:]*d_t/d_x)*(c_w[-1,:]-c_w[-2,:])

    future_wavefield[0,:] = future_wavefield[1,:] #offenes Ende oben

    return future_wavefield

def wellen_simulieren_ohne_wände(u, v, n_t, d_t, d_x, d_y):
    if numpy.max(v)*d_t/d_x > 1: #Sicherstellen, dass die Simulation möglich ist
        print(numpy.max(v)*d_t/d_x)
        print("v*d_t/d_x darf nicht größer als 1 sein. Sonst ist v schneller, als die Elemente nacheinander berechnet werden können.")
        
    for i in range(n_t-2):
        print(str(i+2)+" von "+str(n_t), end='\r')
        u[:,:, i+2] = future_wavefield_ohne_wände(v, u[:, :, i], u[:, :, i+1], d_t, d_x, d_y)
    return u

def future_wavefield_ohne_wände(v, u_p, c_w, d_t, d_x, d_y):
    u_c = numpy.pad(c_w, ((1, 1), (1, 1)), mode='constant')  # Zero-padding

    x_diff = ((v * d_t / d_x) ** 2) * (u_c[2:, 1:-1] - 2 * u_c[1:-1, 1:-1] + u_c[:-2, 1:-1])
    y_diff = ((v * d_t / d_y) ** 2) * (u_c[1:-1, 2:] - 2 * u_c[1:-1, 1:-1] + u_c[1:-1, :-2])

    future_wavefield = 2 * u_c[1:-1, 1:-1] - u_p + x_diff + y_diff
    
    future_wavefield[:,-1] = c_w[:,-1] - (v[:,-2]*d_t/d_x)*(c_w[:,-1]-c_w[:,-2])
    future_wavefield[0,:] = c_w[0,:] - (v[1,:]*d_t/d_x)*(c_w[0,:]-c_w[1,:])
    future_wavefield[:,0] = c_w[:,0] - (v[:,1]*d_t/d_x)*(c_w[:,0]-c_w[:,1])
    future_wavefield[-1,:] = c_w[-1,:] - (v[-2,:]*d_t/d_x)*(c_w[-1,:]-c_w[-2,:])

    return future_wavefield

def point_source(amplitude_gauss, mu_x, mu_y, sigma, n_x, n_y):
    seitenverhältnis = n_x/n_y
    x, y = numpy.meshgrid(numpy.linspace(0, 1, n_x), numpy.linspace(0, 1, n_y) )
    
    gauss = amplitude_gauss * numpy.exp(-0.5 * ((x - mu_x) / (sigma))**2 - 0.5 * ((y - mu_y) / (sigma*seitenverhältnis))**2)
    return gauss