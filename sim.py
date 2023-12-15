import numpy as np # (pip install numpy)
import math
import tomllib # Module qui permet d'utiliser les fichiers toml avec plus de facilité
from loguru import logger # Module qui permet d'avoir un beau logging dans la console (pip install loguru)
import matplotlib.pyplot as plt # (pip install matplotlib)

# Toutes les constantes qui peuvent être modifiées sont dans le fichier data.toml
# pour une raison de simplicité pour ne pas devoir revenir dans le code à chaque fois

#======== CONSTANTES ========
FICHIER_TOML = "data.toml"
DENSITE_EAU = 997 #[kg/m³]
#============================

try: #Lecture du fichier toml avec les configurations de la grue
    tomlFile = open(FICHIER_TOML, "rb")
    data = tomllib.load(tomlFile)
    logger.debug("La lecture de data.toml a été effectuée avec succès")
except tomllib.TOMLDecodeError as e:
    logger.error(f"Erreur dans data.toml: {e}")
    exit()

def angles_immersion_soulevement():
    """
    Calcule les angles d'immersion et de soulèvement
    Pré: /
    Post: retourne une liste avec les deux angles
    """
    theta_max_1 = math.atan((data["Barge"]["Hauteur"]-enfoncement())/(data["Barge"]["Longueur"]/2))
    theta_max_2 = math.atan(enfoncement()/(data["Barge"]["Longueur"]/2))
    return min(theta_max_1, theta_max_2) # Retourne l'angle le plus petit qui provoquerait une immersion

def enfoncement():
    """ 
    Calcule la hauteur de l'enfoncement de la plateforme
    Pré:  /
    Post: retourne un int (hauteur de l'enfoncement)
    """
    return masse_totale() / (DENSITE_EAU * data["Barge"]["Longueur"]**2) # Applique la formule de notre modèle de physique

def masse_articulations():
    """
    Calcule la masse totale de toute la grue et éventuellement d'un poids
    pré: ne prend aucun paramètre
    post: retourne la masse totale de la grue 
    """
    masse_articulations_totale = 0 # Initialise la masse
    for i in range(len(data["Articulations"])): # Additionne la masse de chaque articulation
        masse_articulations_totale += data["Articulations"][str(i)]["masse"]
    masse_articulations_totale += data["ContrePoids"]["Masse"] # Ajoute également le contrepoids
    return masse_articulations_totale

def masse_totale():
    """
    Calcule la masse totale de toute la structure
    pré: ne prend aucun paramètre
    post: retourne la masse totale de la structure 
    """
    return masse_articulations() + data["Barge"]["Masse"] # La seule masse à ajouter est la masse de la barge

# =================================================================
# ======================= CENTRES DE MASSES =======================
# =================================================================

def centre_masse_barge():
    """
    Cette fonction calcule le centre de masse dans un plan 2d de la barge.
		(L'axe x est au niveau de l'eau et y milieu de la barge) 
    Pré: /
    Post: Retourne un tuple avec les coordonnées du centre de masse
    """
    logger.info("Début du calcul du centre de masse de la barge")
    CMx = (data["Barge"]["Longueur"])/2
    CMy = (data["Barge"]["Hauteur"])/2
    
    CMx = CMx-data["Barge"]["Longueur"]/2 # Changement de repère pour le mettre au centre de la barge et à hauteur de l'enfoncement
    CMy = CMy-enfoncement()
    
    CM = (CMx, CMy)

    logger.success(f"Le centre de masse de la barge trouvé est {CM}[m]")
    return CM

def centre_masse_grue(angles_articulations):
    """
    Cette fonction calcule le centre de masse dans un plan 2d de la grue.
		(L'axe x est au niveau de l'eau et y milieu de la barge) 
    Pré: Une liste des angles des articulations
    Post: Retourne un tuple avec les coordonnées du centre de masse
    """
    try:
        logger.info("Début du calcul du centre de masse de la grue")
        extr_socle = {"x": data["Articulations"]["0"]["largeur"]/2, "y": data["Articulations"]["0"]["longueur"]} # La longueur du socle est la hauteur (y) et la largeur/2 le centre du socle (x)
        
        Ex = [extr_socle["x"]] # Initialise les coordonnées des extrémités (avec le socle)
        Ey = [extr_socle["y"]]
        
        for i in range(1,len(data["Articulations"])): # Calcule le x,y de l'extrémité de chaque articulations et le met dans Ex/Ey
            try:
                Ex.append(data["Articulations"][str(i)]["longueur"]*math.cos(angles_articulations[i-1]))
                Ey.append(data["Articulations"][str(i)]["longueur"]*math.sin(angles_articulations[i-1])) 
            except IndexError:
                logger.error(f"Le nombre d'angle n'est pas correct, il faut un total de {len(data['Articulations'])-1} angle(s).")
                exit()
        
        Cx = [extr_socle["x"]] # Dans ces listes seront stockées les centres de masse, on y met déjà le centre de masse du socle
        Cy = [extr_socle["y"]/2]
        
        for i in range(1,len(Ex)): # Pour toutes les articulations on calcule le centre de masse (Pas 0 (socle) car déjà calculé dans les tableaux)
            temp_Cx = 0
            temp_Cy = 0

            for u in range(i+1): # Boucle for pour faire la somme et calculer toutes les articulations
                if i == u: # Si on arrive à l'articulation où l'on cherche le CM, on ajoute les (x,y) de toutes les extrémités précédentes + 1/2 de celle qu'on cherche le CM
                    temp_Cx += Ex[u]/2
                    temp_Cy += Ey[u]/2
                    Cx.append(temp_Cx)
                    Cy.append(temp_Cy)
                    break
                temp_Cx += Ex[u] # Sinon on additionne simplement les extrémités
                temp_Cy += Ey[u]

        pos_massX = 0
        pos_massY = 0

        for i in range(len(Cx)): # Boucle qui va appliquer masse * position (numérateur formule CM) sans le contre-poids
            pos_massX += data["Articulations"][str(i)]["masse"] * Cx[i]
            pos_massY += data["Articulations"][str(i)]["masse"] * Cy[i]

        CMx_contre_poids = -data["ContrePoids"]["Longueur"]*math.cos(angles_articulations[0])/2 # Le centre de masse (en x) du contre-poids est moins sa longueur * sin de l'angle de la première articulation, le tout divisé par 2
        CMy_contre_poids = data["Articulations"]["0"]["longueur"] - data["ContrePoids"]["Longueur"]*math.sin(angles_articulations[0])/2 # Le centre de masse (en y) du contre-poids est la hauteur du socle moins la longueur du contre-poids/2

        pos_massX += CMx_contre_poids * data["ContrePoids"]["Masse"]
        pos_massY += CMy_contre_poids * data["ContrePoids"]["Masse"]

        CMx = pos_massX/masse_articulations() # Applique la formule du centre de masse
        CMy = pos_massY/masse_articulations()

        CMx = CMx+data["Barge"]["Declage_grue"]-data["Barge"]["Longueur"]/2 # Changement de repère pour le mettre au centre de la barge et à hauteur de l'enfoncement
        CMy = CMy+data["Barge"]["Hauteur"]-enfoncement()

        CM = (CMx,CMy)
        logger.success(f"Le centre de masse de la grue trouvé est {CM}[m]")
        return CM
    except KeyError as e:
        logger.error(f"Vous avez fait une erreur dans le fichier toml, Python le renseigne l'index {e}")
        exit()

def centre_masse_total(angles_articulations):
    """
    Cette fonction calcule le centre de masse de l'ensemble de la structure
    Pré: Récupère une liste d'angles exprimés en radian contenant n-1 éléments du dictionnaire
    Post: Retourne un tuple avec les coordonnées du centre de masse
    """
    CMgrue = centre_masse_grue(angles_articulations)
    CMbarge = centre_masse_barge()
    
    CMx = (CMbarge[0] * data["Barge"]["Masse"] + CMgrue[0] * masse_articulations())/masse_totale() # Applique la formule générale du centre de masse pour les composantes x et y
    CMy = (CMbarge[1] * data["Barge"]["Masse"] + CMgrue[1] * masse_articulations())/masse_totale()
    
    CM = (round(CMx, 4), round(CMy, 4)) # Arrondis à 4 chiffres après la virgule
    logger.success(f"Le centre de masse total est {CM}[m]")
    return CM

# =================================================================

def base_trapezes(theta):
    """
    Calcule la base des trapèzes du volume immergé
    Pré: Prend l'angle théta
    Post: Retourne les deux côtés du trapèze
    """
    hr = enfoncement() + data["Barge"]["Longueur"]/2 * math.tan(theta)
    hl = enfoncement() - data["Barge"]["Longueur"]/2 * math.tan(theta)
    return (hr, hl)

def centre_poussee(theta):
    """
    Calcule la position du centre de poussée de la barge en fonction de théta
    Pré: L'angle d'inclinaison en radian de la barge (float)
    Post: Retourne position en x du centre de poussée par rapport à l'origine du repère.
    """
    hr = base_trapezes(theta)[0]
    hl = base_trapezes(theta)[1]
    L = data["Barge"]["Longueur"]
    x_sec = (L * (2*hl + hr))/(3 * (hl + hr)) - L/2     # -L/v et -enfoncement() font le changement de repère
    y_sec = (hr**2 + hr * hl + hl**2)/(3 * (hl + hr)) - enfoncement()
    x_prim = x_sec * math.cos(theta) - y_sec * math.sin(theta) # Effectue la rotation pour l'axe x
    y_prim = x_sec * math.sin(theta) + y_sec * math.cos(theta) # Effectue la rotation pour l'axe y
    return (x_prim, y_prim)

def inertie():
    """
    Calcule l'inertie de la barge
    Pré: /
    Post: Retourne un float qui donne l'inertie de la barge
    """
    return (masse_totale()*(data["Barge"]["Longueur"])**2+data["Barge"]["Hauteur"]**2)/12

def convertisseur_tracker(fichier):
    """
    Fonction qui retourne les angles théta en fonction du temps du fichier tracker
    Pré: Prend le nom du fichier
    Post: Retourne le temps avec son angle
    """
    read_data = np.loadtxt(fichier) # Charge les données du fichier
    (t,x1,y1,x2,y2,x3,y3) = read_data.T # Dézippe les colonnes du tableau de données dans des variables distinctes (t, x1, y1, x2, y2, x3, y3)
    len_time = len(t) # Obtient la longueur du vecteur temps

    vec1 = np.zeros((len_time,2))   # Initialise deux vecteurs pour stocker les composantes des vecteurs entre les points
    vec2 = np.zeros((len_time,2))
    cos_angle = np.zeros_like(t)    # Initialise un tableau pour stocker les cosinus des angles entre les vecteurs
    theta = np.zeros_like(t)        # Initialise un tableau pour stocker les angles theta
    
    for i in range(len_time): # Boucle à travers les points pour calculer les vecteurs, les cosinus et les angles theta
        vec1[i,:] = [x2[i]-x1[i],y2[i]-y1[i]]
        vec2[i,:] = [x3[i]-x1[i],y3[i]-y1[i]]
        cos_angle[i] = np.dot(vec1[i,:],vec2[i,:])/( np.linalg.norm(vec1[i,:]) * np.linalg.norm(vec2[i,:]))
        theta[i] = math.acos(cos_angle[i])*(180/math.pi)-90 
    return t, theta # Retourne le temps et les angles theta calculés

step = 0.0001 # Le t+1 infinitésimal calculé
end = 5.5 # Le t maximal représenté sur le graphe
t = np.arange(0, end, step) # Le tableau numpy avec tout les tmps 
theta = np.empty_like(t)
omega = np.empty_like(t)
accelerationAngulaire = np.empty_like(t)
max_angle_array = np.full_like(t, -angles_immersion_soulevement())

Eg = np.empty_like(t) # Variables pour le graphique des énergies
Ec = np.empty_like(t)
Ek = np.empty_like(t)
Ea = np.empty_like(t)
Etot = np.empty_like(t)

CMgrue = centre_masse_grue(data["Angles"]) # Calcule le centre de masse de la grue en fonction des angles d'inclinaison des articulations
CMtotal = centre_masse_total(data["Angles"]) # Calcule le centre de masse total en fonction des angles d'inclinaison des articulations

theta[0] = 0 # S'initialise en theta
dt = step

logger.info(f"La masse totale de la grue est de {masse_totale()}") # Affiche les données importantes
logger.info(f"L'inertie de la grue est de {inertie()}")
logger.info(f"La constante d'amortissement est {data['Barge']['ConstanteAmortissement']}")

logger.debug("Début de la méthode d'Euler explicite")
for x in range(len(t)-1): # Boucle qui permet pour chaque dt de calculer les thétas
    dt = step
    
    centrePousseeX = centre_poussee(theta[x])[0] 
    centrePousseeY = math.cos(theta[x]) * centre_poussee(theta[x])[1]
    centreGraviteX = math.sin(theta[x]) * CMtotal[1] # Prend en compte la rotation du centre de masse total
    centreGraviteY = math.cos(theta[x]) * CMtotal[1]

    coupleRedressement = masse_totale() * 9.81 * abs(centrePousseeX - centreGraviteX) # Cr
    coupleDestabilisateur = -masse_articulations() * 9.81 * CMgrue[0] # Ca
    totalCouples = coupleDestabilisateur + coupleRedressement

    accelerationAngulaire[x+1] = (-data["Barge"]["ConstanteAmortissement"] * omega[x] + totalCouples) / inertie() # Application de la méthode d'Euler explicite
    omega[x+1] = omega[x] + accelerationAngulaire[x+1] * dt
    theta[x+1] = theta[x] + omega[x+1] * dt

    Eg[x+1] = masse_totale() * 9.81 * (centreGraviteY-CMtotal[1])               # E Flotteur
    Ec[x+1] = -masse_totale() * 9.81 * (centrePousseeY-centre_poussee(0)[1])    # E Poussée
    Ek[x+1] = enfoncement() * omega[x]**2/2                                     # E Cinétique
    Ea[x+1] = -coupleDestabilisateur * theta[x]                                 # E Charge
    Etot[x+1] = Eg[x] + Ec[x] + Ek[x] + Ea[x]                                   # E Totale

logger.debug("Fin de la méthode d'Euler") 

tracker = convertisseur_tracker("data.txt")

plt.figure(figsize=(9, 7)) # Donne une taille plus grande à la fenêtre du graphique
plt.subplot(3, 1, 1) # Premier graphique
plt.plot(t, np.degrees(theta), label="θ", color="green", linewidth=1)
plt.plot(tracker[0], tracker[1]+67, label="tracker", color="pink")
plt.plot(t, np.degrees(max_angle_array), "--", label="θ max", color="red", linewidth=1)
plt.xlabel("temps (s)")
plt.ylabel("angle (°)")

plt.title("Angle/temps")
plt.annotate(round(theta[-1], 3), (27, theta[-1] + 0.7))
plt.legend(prop={'size': 6}, loc=4)

plt.subplot(3, 1, 2) # Deuxième graphique
plt.plot(t, np.rad2deg(omega), label="ω", color="green", linewidth=1)
plt.xlabel("temps (s)")
plt.ylabel("vitesse (°/s)")

plt.title("Vitesse/temps")
plt.legend(prop={'size': 6})

plt.subplot(3, 1, 3) # Troisième graphique
plt.plot(t, np.rad2deg(accelerationAngulaire), label="α", color="green", linewidth=1)
plt.xlabel("temps (s)")
plt.ylabel("accélération (°/s²)")

plt.title("Accélération/temps")
plt.legend(prop={'size': 6})

plt.tight_layout() # Ajuste l'espacement entre les graphiques pour éviter un chevauchement
plt.show() # Affiche les trois graphiques

plt.figure(5)
plt.plot(np.rad2deg(omega), np.rad2deg(theta), color="green")
plt.xlabel("theta (°)")
plt.ylabel("vitesse angulaire (°/s)")
plt.title("Diagramme de phase")
plt.grid(True) # Affiche un cadrillage
plt.show() # Affiche le diagramme de phase

plt.figure(2)
plt.plot(t, Eg, label="E Flotteur", color="#06668C")
plt.plot(t, Ec, label="E Poussée", color="#E9C46A")
plt.plot(t, Ek, label="E Cinétique", color="#A4BD01")
plt.plot(t, Ea, label="E Charge", color="#BD3100")
plt.plot(t, Etot, label="E Totale", color="#80586D")
plt.xlabel("temps (s)")
plt.ylabel("énergie (J)")
plt.title("Energie/temps")
plt.legend(prop={'size': 7},loc='upper right')
plt.show() # Affiche le graphique des énergies