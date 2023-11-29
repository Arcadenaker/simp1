import numpy as np
import math
import tomllib
from loguru import logger
import matplotlib.pyplot as plt

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

def gravite(masse):
    """
    Calcule la force de garvité de la masse donnée
    Pré: prend la masse en paramètre (int)
    Post: retourne la force en [N] (int)
    """
    return masse * 9.81

def angles_immersion_soulevement():
    """
    Calcule les angles d'immersion et de soulèvement
    Pré: /
    Post: retourne une liste avec les deux angles
    """
    theta = [math.atan((data["Barge"]["Hauteur"]-enfoncement())/(data["Barge"]["Longueur"]/2)), math.atan(enfoncement()/(data["Barge"]["Longueur"]/2))]
    return min(theta)

def enfoncement():
    """ 
    Calcule la hauteur de l'enfoncement de la plateforme
    Pré:  /
    Post: retourne un int (hauteur de l'enfoncement)
    """
    return masse_totale() / (DENSITE_EAU * data["Barge"]["Longueur"]**2)

def masse_totale():
    """
    Calcule la masse totale de toute la structure
    pré: ne prend aucun paramètre
    post: retourne la masse totale de la structure 
    """
    sum = 0
    for m in range(1,len(data["Articulations"])-1):
        sum += int(data["Articulations"][str(m)]["masse"])
    sum += data["Barge"]["Masse"] + data["ContrePoids"]["Masse"]
    return sum

def masse_articulations():
    """
    Calcule la masse totale de toute la grue et éventuellement d'un poids
    pré: ne prend aucun paramètre
    post: retourne la masse totale de la grue 
    """
    masse_articulations_totale = 0
    for i in range(len(data["Articulations"])):
        masse_articulations_totale += data["Articulations"][str(i)]["masse"]
    masse_articulations_totale += data["ContrePoids"]["Masse"]
    return masse_articulations_totale

# ======================= CENTRES DE MASSES =======================

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

def centre_masse_grue(angles):
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
                Ex.append(data["Articulations"][str(i)]["longueur"]*math.cos(angles[i-1]))
                Ey.append(data["Articulations"][str(i)]["longueur"]*math.sin(angles[i-1])) 
            except IndexError:
                logger.error(f"Le nombre d'angle n'est pas correct, il faut un total de {len(data['Grue'])-1} angle(s).")
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

        CMx_contre_poids = -data["ContrePoids"]["Longueur"]*math.cos(angles[0])/2 # Le centre de masse (en x) du contre-poids est moins sa longueur * sin de l'angle de la première articulation, le tout divisé par 2
        CMy_contre_poids = data["Articulations"]["0"]["longueur"] - data["ContrePoids"]["Longueur"]*math.sin(angles[0])/2 # Le centre de masse (en y) du contre-poids est la hauteur du socle moins la longueur du contre-poids/2

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

def centre_poussee(theta):
    """
    Calcule la position du centre de poussée de la barge en fonction de théta
    Pré: L'angle d'inclinaison en radian de la barge (float)
    Post: Retourne position en x du centre de poussée par rapport à l'origine du repère.
    """
    return (math.sin((theta)) * data["Barge"]["Longueur"]**2) / (12 * enfoncement() * ((math.cos(theta)) ** 2))

def inertie():
    """
    Calcule l'inertie de la barge
    Pré: /
    Post: Retourne un float qui donne l'inertie de la barge
    """
    return (masse_totale()*(data["Barge"]["Longueur"])**2+data["Barge"]["Hauteur"]**2)/12

step = 0.0001 # Le t+1 infinitésimal calculé
end = 5 # Le t maximal représenté sur le graphe
t = np.arange(0, end, step)
theta = np.empty_like(t)
omega = np.empty_like(t)
accelerationAngulaire = np.empty_like(t)
max_angle_array = np.full_like(t, -angles_immersion_soulevement())

CMgrue = centre_masse_grue(data["Angles"])
CMtotal = centre_masse_total(data["Angles"])

theta[0] = 0 # S'initialise en theta
dt = step

for x in range(len(t)-1):
    dt = step
    centrePoussee = centre_poussee(theta[x])
    centreGravite = math.sin(theta[x]) * CMtotal[0]

    coupleRedressement = gravite(masse_totale()) * abs(centrePoussee - centreGravite)
    coupleDestabilisateur = -gravite(masse_articulations()) * CMgrue[0]
    totalCouples = coupleDestabilisateur + coupleRedressement

    accelerationAngulaire[x+1] = (-data["Barge"]["ConstanteAmortissement"] * omega[x] + totalCouples) / inertie() # Application de la méthode de Gauss
    omega[x+1] = omega[x] + accelerationAngulaire[x+1] * dt
    theta[x+1] = theta[x] + omega[x+1] * dt


# Premier graphique
plt.subplot(3, 1, 1)
plt.plot(t, np.degrees(theta), label="θ", color="green", linewidth=1)
plt.plot(t, np.degrees(max_angle_array), "--", label="θ min", color="purple", linewidth=1)
plt.xlabel("temps (s)")
plt.ylabel("angle (°)")
plt.title("Angle/temps")
plt.annotate(round(theta[-1], 3), (27, theta[-1] + 0.7))
plt.legend(prop={'size': 6})

# Deuxième graphique
plt.subplot(3, 1, 2)
plt.plot(t, np.rad2deg(theta), label="ω", color="green", linewidth=1)
plt.xlabel("temps (s)")
plt.ylabel("vitesse (°/s)")
plt.title("Vitesse/temps")
plt.legend(prop={'size': 6})

# Troisième graphique
plt.subplot(3, 1, 3)
plt.plot(t, np.rad2deg(accelerationAngulaire), label="α", color="green", linewidth=1)
plt.xlabel("temps (s)")
plt.ylabel("accélération (°/s^2)")
plt.title("Accélération/temps")
plt.legend(prop={'size': 6})

# Ajuster l'espacement entre les graphiques pour éviter le chevauchement
plt.tight_layout()

# Afficher les graphiques
plt.show()