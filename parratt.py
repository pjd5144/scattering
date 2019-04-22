def parratt_reflectivity(wavelength,n_substrate,n_film,t_film,alpha_i):
    k0 = 2*np.pi/wavelength
    num_slice = 150.
    slice_thick = t_film/num_slice
    layers = np.arange(0,-t_film-slice_thick,-slice_thick)
    n_layers = np.zeros((len(layers),1),dtype=np.complex_)
    X = np.zeros((len(layers),1),dtype=np.complex_)

    # layer assignments
    n_layers[0] = 1 + 0.j
    n_layers[1:-1] = n_film
    n_layers[-1] = n_substrate

    Rf = np.zeros((len(alphai),1),dtype=np.complex_)
    idx = 0
    for k in alphai:
        # z-component of wavevector
        kz = k0*np.sqrt(n_layers**2-np.cos(k)**2)
        r = (kz[0:-1] - kz[1:len(n_layers)+1])/(kz[0:-1] + kz[1:len(n_layers)+1])

        for i in range(len(n_layers)-2,-1,-1):
            X[i] = (np.exp(-2.j*kz[i]*layers[i]) *
                    (r[i]+X[i+1]*np.exp(2.j*kz[i+1]*layers[i])) /
                    (1+r[i]*X[i+1]*np.exp(2.j*kz[i+1]*layers[i])))
        Rf[idx] = X[0]
        idx += 1
    return Rf
