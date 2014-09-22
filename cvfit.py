def a(n):
	if n <= 300:
		nbdec = 5.0
	elif n <= 1000:
		nbdec = 8.0
	elif n <= 3000:
		nbdec = 10.0
	else:
		nbdec = 12.0
	return nbdec

def b(n):
	return 5.0 + (n > 300) * 3 + (n > 1000) * 2 + (n > 3000) * 2
