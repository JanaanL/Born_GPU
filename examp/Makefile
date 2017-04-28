

data.json: data.P
	Data json=data.P	


wavelet.json: wavelet.P
	Wavelet json=wavelet.P


vel.json: vel.P
	Vel json=vel.P

mig.json: vel.json wavelet.json data.json mig.P
	RTM3D json=mig.P 
	
