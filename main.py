###########
## Dados ##
###########

# prob de contaminação/comutação entre estados
d = 5
# prob de diaria de recuperação r = 1-d
r = 1
# risco de transmição em area com alta densidade de mosquitos - h = [0.12, 11.2]
h = 0.5
# risco de transmição em area com baixa densidade de mosquitos - l = [0.09, 1.02]
l = 0.5
# dias infectado
D = 5
# taxa de picadas de mosquito - a = [0.6, 1.2]
al = 0.6
ah = 1.2
# infecção mosquito -> humano - b = [0.5, 0.9]
# infecção humano -> mosquito - c = [0.5, 0.9]
bl = cl = 0.5
bh = ch = 0.9
# periodo de incubação t = [7, 25]
tl = 13
#th = 25
# expectativa de vida do mosquito/taxa de morte - 1/u = [0.83, 0.96]
ul = 0.09
uh = 0.09
# densidade de mosquitos em 2010
M = 2.24
# densidade de humandos em 2010
H = 1215.69

# numero de Euler
e = 2.718

# função para recuperar tabela do txt
# assume que é uma matriz quadrada
def recuperaTabela():

    tabela = []
    arquivo = open('tabela.txt', 'r') # abre arquivo

    while True:# enquanto !EOF
        temp = arquivo.readline() # le linha do txt
        if(temp == ''):
            break
        else:
            vetorTemp = temp.split()

            # percorre o vetor transformando os itens em float
            for i in range(len(vetorTemp)):
                vetorTemp[i] = float(vetorTemp[i])

            # adiciona itens na tabela
            tabela.append(vetorTemp)

    # fecha arquivo
    arquivo.close()
    return tabela

# estimando yh/yl
def casosSecundarios(a,b,c,u,t):
    y = a*(M/H)*a*((e**(-u*t))/u)*b*c
    return y

# Estimando Dh/Dl
def tempoPassado(i):
    Di = ((d+i)/D)/(1/(D*(2*d + 1/D)))
    return Di

# chapman-kolmogorov
def Chapman_Kolmogorov(tabela,step, i, j,tam): # funcao principal
    probVisit = 0.0
    if step == 1:
        return tabela[i][j]
    else:
        for k in range(tam):
            probVisit += tabela[i][k] * Chapman_Kolmogorov(tabela,step-1,k,j,tam)
    return probVisit

######################
## rotina principal ##
######################

tabela=recuperaTabela()
tam = len(tabela) # ordem da tabela e tamanho do vetor (numero de colunas)
# exibe tabela
print("tabela: ")
for i in range(tam):
    print()
    for j in range(tam):
        print(tabela[i][j],end=' | ')
print("\n")

print("Casos secundarios em area com alta densidade de vetores:")
print(casosSecundarios(ah,bh,ch,uh,tl))
print()
print("Casos secundarios em area com baixa densidade de vetores:")
print(casosSecundarios(al,bl,cl,ul,tl))
print()

print("Tempo gasto em area com alta densidade de vetores")
print(tempoPassado(h))
print()
print("Tempo gasto em area com baixa densidade de vetores")
print(tempoPassado(l))
print()

print("Potencial perigo de epidemia")
print((tempoPassado(h)*casosSecundarios(ah,bh,ch,uh,tl)) + (tempoPassado(l)*casosSecundarios(al,bl,cl,ul,tl)))

## elementos estacionarios ##
print()
A = (d+r)/(3*(2*d+r)*(e*r+d))
piDelta = r/(2*r + d)
piYh = A*(0.3*(d-r) + r + 2*d + 0.7*(2*r+d))
piYl = A*(0.3*(2*r-d) + r + 2*d + 0.7*(d-r))

print("Elementos de distribuição estacionaria")
print("piDelta:",piDelta)
print("pi Yh:",piYh)
print("pi Yl:",piYl)


print("\nmatriz depois de 6 anos (2016)")
for i in range(tam):
    print()
    for j in range(tam):
        result = Chapman_Kolmogorov(tabela,6,i,j,tam)
        print(result,end=' | ')