import chemev

chab = chemev.imf.Chabrier()
print 'Number of Stars with M>#:'
print '0:',chab.cum_num(0)
print '1:',chab.cum_num(1)
print '2:',chab.cum_num(2)
print '3:',chab.cum_num(3)
print '10:',chab.cum_num(10)
print '0,1,2,3,10:',chab.cum_num([0,1,2,3,10])

print 'Stellar Masses:'
print '4:',chab.cum_mass(4)
print '7:',chab.cum_mass(7)
print '9:',chab.cum_mass(9)
print '0:',chab.cum_mass(0)
print '4,7,9,0:',chab.cum_mass([4,7,9,0])

