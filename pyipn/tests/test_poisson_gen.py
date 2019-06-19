from pyipn.poisson_gen import source_poisson_generator, background_poisson_generator

def test_source_poisson_gen():
	d1 = source_poisson_generator.__wrapped__(0.0, 1, 1, 0.5, 0.5, 0.5)
	d2 = source_poisson_generator.__wrapped__(-50, 145.5, 100, 0, 0.5, 10)
	d2 = source_poisson_generator.__wrapped__(-50, 145.5, 0, 0, 0.5, 10)

def test_background_poisson_gen():
	b1 = background_poisson_generator.__wrapped__(0.0, 1, 1, 0.5)
	b2 = background_poisson_generator.__wrapped__(-50, 145.5, 100, 0) 
