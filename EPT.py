def EPT(j_me):
    ept = open("EarthPeriodicTerms.txt", 'r')
    categories = ept.readline()

    def splitter(num):
        l = []
        for i in range(num):
            line = ept.readline()
            line = line.split(' ')
            length = len(line)
            acc = 0
            #removes the four '' in the list
            while acc<length:
                if line[acc] == '':
                    del line[acc]
                    length -= 1
                else:
                    break
            #removes /n in list
            last = list(line[-1])
            last = last[:-1]
            last = ''.join(last)
            line[-1] = last
            
            if len(line) == 5:
                l.append(line[1:])
            else:
                l.append(line)
        return l
    l0_t = splitter(64)
    l1_t = splitter(34)
    l2_t = splitter(20)
    l3_t = splitter(7)
    l4_t = splitter(3)
    l5_t = splitter(1)
    def eq(l, j): #j = jme
        from math import cos
        ln = 0
        for i in l:
            a, b, c = float(i[1]), float(i[2]), float(i[3])
            l_i = a * cos(b+(c*float(j)))
            ln += l_i
        return ln
    l0 = eq(l0_t, j_me)
    l1 = eq(l1_t, j_me)
    l2 = eq(l2_t, j_me)
    l3 = eq(l3_t, j_me)
    l4 = eq(l4_t, j_me)
    l5 = eq(l5_t, j_me)   
    return [l0, l1, l2, l3, l4, l5]

############################################################

def B_EPT(j_me):
    b_ept = open("B_EarthPTerms.txt", 'r')
   
    def splitter(num):
        b = []
        for i in range(num):
            line = b_ept.readline()
            line = line.split(' ')
            length = len(line)
            acc = 0
            #removes the four '' in the list
            while acc<length:
                if line[acc] == '':
                    del line[acc]
                    length -= 1
                else:
                    break
            #removes /n in list
            last = list(line[-1])
            last = last[:-1]
            last = ''.join(last)
            line[-1] = last
            
            if len(line) == 5:
                b.append(line[1:])
            else:
                b.append(line)
        return b
    b0_t = splitter(5)
    b1_t = splitter(2)

    def eq(B, j): #j = jme
        from math import cos
        Bn = 0
        for i in B:
            a, b, c = float(i[1]), float(i[2]), float(i[3])
            B_i = a * cos(b+(c*float(j)))
            Bn += B_i
        return Bn
    b0 = eq(b0_t, j_me)
    b1 = eq(b1_t, j_me)
    return [b0, b1]

############################################################

def R_EPT(j_me):
    r_ept = open("R_EarthPTerms.txt", 'r')
   
    def splitter(num):
        r = []
        for i in range(num):
            line = r_ept.readline()
            line = line.split(' ')
            length = len(line)
            acc = 0
            #removes the four '' in the list
            while acc<length:
                if line[acc] == '':
                    del line[acc]
                    length -= 1
                else:
                    break
            #removes /n in list
            last = list(line[-1])
            last = last[:-1]
            last = ''.join(last)
            line[-1] = last
            
            if len(line) == 5:
                r.append(line[1:])
            else:
                r.append(line)
        return r
    r0_t = splitter(40)
    r1_t = splitter(10)
    r2_t = splitter(6)
    r3_t = splitter(2)
    r4_t = splitter(1)
   
    def eq(R, j): #j = jme
        from math import cos
        Rn = 0
        for i in R:
            a, b, c = float(i[1]), float(i[2]), float(i[3])
            R_i = a * cos(b+(c*float(j)))
            Rn += R_i
        return Rn
    r0 = eq(r0_t, j_me)
    r1 = eq(r1_t, j_me)
    r2 = eq(r2_t, j_me)
    r3 = eq(r3_t, j_me)
    r4 = eq(r4_t, j_me)
    return [r0, r1, r2, r3, r4]

############################################################

def nutation(j_ce, X0, X1, X2, X3, X4):
    from math import sin, cos, radians
    nu = open("NutationTerms.txt", 'r')
    categories = nu.readline()  
    n = []
    for i in range(63):
        line = nu.readline()
        line = line.split(' ')
        last = list(line[-1])
        last = last[:-1]
        last = ''.join(last)
        line[-1] = last
        n.append(line)
    psi_sum = 0
    epsilon_sum = 0
    for i in n:
        y0, y1, y2, y3, y4 = float(i[0]), float(i[1]), float(i[2]), float(i[3]), float(i[4])
        a, b, c, d = float(i[5]), float(i[6]), float(i[7]), float(i[8])
        sumxy = (float(X0)*y0)+(float(X1)*y1)+(float(X2)*y2)+(float(X3)*y3)+(float(X4)*y4)
        psi_i = (a+b*float(j_ce))*sin(radians(sumxy))
        epsilon_i = (c+d*float(j_ce))*cos(radians(sumxy))
        psi_sum += psi_i
        epsilon_sum += epsilon_i
    psi_d = psi_sum/36000000
    epsilon_d = epsilon_sum/36000000
    return [psi_d, epsilon_d]
    
        
