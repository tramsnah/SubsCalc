import numpy as np

def vanopstal(x,y,c,k,xpm,ypm,comp):
    '''
    Functie Van Opstal (1974) The effect of base-rock rigidity on subsidence 
    due to reservoir compaction. Proc. 3rd Congr. Int. Soc. Rock Mech., 
    Denver Colorado, sept. 1974, vol 2, part B.
    
    geschreven door Karin van Thienen-Visser, TNO-AGE
    
    Te gebruiken: subsidence=vanopstal(x,y,c,k,xpm,ypm,comp)
    
    Input:
    xpm, ypm zijn de coordinaten van de punten waarop de bodemdaling uitgerekend wordt
            (numpy arrays)
    x and y zijn de coordinaten van de gridcellen in het reservoir
            (numpy arrays)
    c is de diepte van het reservoir. Dit is een numpy array met dezelfde lengte als x,
            maar c kan ook een getal zijn
    k is de diepte van het rigide basement
    comp is de compactie in elk van de gridcellen van het reservoir,
            vergelijkbaar met x,y,c, input
            De eenheid van 'comp' is m3, dus typisch;
                comp=dx*dy*h*dp*Cm
    Output:
    numpy array met zelfde afmetingen als xpm. Bodemdaling is positief geformuleerd.
    '''

    ##start van de functie
    # allocatie van parameters
    subsidence=np.zeros(len(xpm))

    # definieer parameters
    v=0.25      # v=Poisson's ratio, v=0.25 wordt standaard aangenomen in Van Opstal
    a11=1.0
    a12=-0.778
    a21=2.0
    a22=2.80    # Van Opstal (1974) constanten voor v=0.25
    b11=-0.2
    b21=4.0     # Van Opstal (1974) constanten voor v=0.25   

    # bereken van Opstal(1974) termen die alleen afhankelijk zijn van k en c  
    a1=(a11*(a21*k-c)-2*a11*k)
    a2=2*a11*k-(3-4*v)**2*a11*(a21*k+c)
    a3=(3-4*v)*a11*(a21*k+2*k+c)
    a4=(3-4*v)*a11*(a21*k+2*k-c)
    a5=6*a11*k*(a21*k-c)**2
    a6=6*a11*k*(a21*k+c)*(6*k-(a21*k+c))
    a7=60*a11*k**2*(a21*k+c)**3
    a8=a12*(a22*k-c)-2*a12*k
    a9=2*a12*k-(3-4*v)**2*a12*(a22*k+c)
    a10=(3-4*v)*a12*(a22*k+2*k+c)
    aa11=(3-4*v)*a12*(a22*k+2*k-c)
    aa12=6*a12*k*(a22*k-c)**2
    a13=6*a12*k*(a22*k+c)*(6*k-(a22*k+c))
    a14=60*a12*k**2*(a22*k+c)**3

    b1=b11*k
    b2=(3-4*v)**2*b11*k
    b3=(3-4*v)*b11*k
    b4=3*b11*k*(b21*k-c)*(6*k-(b21*k-c))
    b5=3*b11*k*(b21*k+c)*(6*k-(3-4*v)**2*(b21*k+c))-36*b11*k**3
    b6=3*(3-4*v)*b11*k*(b21*k+2*k-c)**2
    b7=3*(3-4*v)*b11*k*(b21*k+2*k+c)**2
    b8=30*b11*k**2*(b21*k-c)**3
    b9=30*b11*k**2*(b21*k+c)**2*(12*k-(b21*k+c))
    b10=420*b11*k**3*(b21*k+c)**4

    ## Berekening invloedsfunctie over alle punten waarover de bodemdaling wordt berekend 
    ## De loop over xpm is expliciet, die over x impliciet
    lx=len(xpm)
    #if(isinstance(c*1.0,float)):
    #    yyy=c
    #else:
    #    yyy=max(c)
    for i in range(lx):
        r=np.sqrt((xpm[i]-x)**2+(ypm[i]-y)**2) # afstand van alle gridcellen in het reservoir t.o.v. de punten aan het oppervlak
        
        # Geertsema oplossing
        u=((1-v)/(2*np.pi))*((2*c)/((r**2+c**2)**(3/2))) 
        #xxx=min(r)
        #if(xxx>4*yyy):   
        #if (max(u)<2e-10):        
        #    subsidence[i]=0
        #    #print(max(u))
        #    continue
            
        # Opstal (1974) termen die afhankelijk zijn van de straal
        n1=(r**2+(a21*k-c)**2)
        n2=(r**2+(a21*k+c)**2)
        terma1=+(a1)/(n1**(3/2))
        terma2=+(a2)/(n2**(3/2))
        terma3=-(a3)/((r**2+(a21*k+2*k+c)**2)**(3/2))
        terma4=+(a4)/((r**2+(a21*k+2*k-c)**2)**(3/2))
        terma5=+(a5)/(n1**(5/2))
        terma6=+(a6)/(n2**(5/2))
        terma7=-(a7)/(n2**(7/2))
        
        n1=(r**2+(a22*k-c)**2) 
        n2=(r**2+(a22*k+c)**2)
        
        terma8=+(a8)/(n1**(3/2))
        terma9=+(a9)/(n2**(3/2))
        terma10=-(a10)/((r**2+(a22*k+2*k+c)**2)**(3/2))
        terma11=+(aa11)/((r**2+(a22*k+2*k-c)**2)**(3/2))
        terma12=+(aa12)/(n1**(5/2))
        terma13=+(a13)/(n2**(5/2))
        terma14=-(a14)/(n2**(7/2))
        
        sum1a=terma1+terma2+terma3+terma4+terma5+terma6+terma7
        sum1b=terma8+terma9+terma10+terma11+terma12+terma13+terma14
        sum1=sum1a+sum1b
        
        n1=(r**2+(b21*k-c)**2)
        n2=(r**2+(b21*k+c)**2)
        n3=(r**2+(b21*k+2*k-c)**2)
        n4=(r**2+(b21*k+2*k+c)**2)
        termb1=-(b1)/(n1**(3/2))
        termb2=+(b2)/(n2**(3/2))
        termb3=-(b3)/(n3**(3/2))
        termb4=+(b3)/(n4**(3/2))
        termb5=-(b4)/(n1**(5/2))
        termb6=+(b5)/(n2**(5/2))
        termb7=+(b6)/(n3**(5/2))
        termb8=-(b7)/(n4**(5/2))
        termb9=+(b8)/(n1**(7/2))
        termb10=+(b9)/(n2**(7/2))
        termb11=-(b10)/(n2**(9/2))
        
        sum2=termb1+termb2+termb3+termb4+termb5+termb6+termb7+termb8+termb9+termb10+termb11
        
        # Van Opstal (1974) correctie op Geertsema
        du=(1-v)/(2*np.pi)*(sum1+sum2)  # bodemdaling is positief geformuleerd 
        
        # Totale nucleus oplossing (Geertsema met Van Opstal correctie)
        g=u+du                      
        
        # Invloedsfunctie op alle punten aan het oppervlak waarover bodemdaling
        # wordt berekend
        subsidence[i]=sum(comp*g)
        
    # einde van loop

    return subsidence
    # einde van Van Opstal functie