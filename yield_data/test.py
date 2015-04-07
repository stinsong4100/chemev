import chemev

chab = chemev.imf.Chabrier()
print chab.cum_num(0)
print chab.cum_num(1)
print chab.cum_num(2)
print chab.cum_num(3)
print chab.cum_num(10)
print chab.cum_num([0,1,2,3,10])

print chab.cum_mass(4)
print chab.cum_mass(7)
print chab.cum_mass(9)
print chab.cum_mass(0)
print chab.cum_mass([4,7,9,0])

