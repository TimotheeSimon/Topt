import os

base_dir = os.path.join(
	os.path.dirname( os.path.realpath( __file__ ) ),
	os.path.join( 'data')
	)

leapseconds_kernel = os.path.join( base_dir, 'latest_leapseconds.tls' )
de432s = os.path.join( base_dir, 'de432s.bsp' )
jup365 = os.path.join( base_dir, 'jup365.bsp' )
