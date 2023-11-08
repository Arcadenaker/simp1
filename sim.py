import numpy as np
import math
import tomllib
from loguru import logger
import matplotlib.pyplot as plt

#======== CONSTANTES ========
FICHIER_TOML = "data.toml"
DENSITE_EAU = 997 #[kg/m³]
#============================

#Lecture du fichier toml avec les configurations de la grue
try:
    tomlFile = open(FICHIER_TOML, "rb")
    data = tomllib.load(tomlFile)
    logger.debug("La lecture de data.toml a été effectuée avec succès")
except tomllib.TOMLDecodeError as e:
    logger.error(f"Erreur dans data.toml: {e}")
    exit()



def centre_masse_grue(articulations, angles):
    """
    Calcule le centre de masse dans un plan 2d de la grue dans un axe oxy dont le centre est le début de la base
    Pré: Récupère un dictionnaire sous la forme {'0': {'longueur': 0, 'largeur': 0, 'masse': 0}...}
            et une liste d'angles exprimés en radian contenant n-1 éléments du dictionnaire
    Post: retourne un tuple décrivant la position x,y du centre de masse
    """
    try:
        logger.info("Début du calcul du centre de masse de la grue")
        extr_socle = {"x": articulations["0"]["largeur"]/2, "y": articulations["0"]["longueur"]/2}
        # Initialise les coordonnées des extrémités (avec le socle)
        Ex = [extr_socle["x"]]
        Ey = [extr_socle["y"]]
        # Calcule le x,y de l'extrémité de chaque articulations et le met dans Ex/Ey
        for i in range(1,len(articulations)):
            try:
                Ex.append(articulations[str(i)]["longueur"]*math.cos(angles[i-1]))
                Ey.append(articulations[str(i)]["longueur"]*math.sin(angles[i-1])) 
            except IndexError:
                logger.error(f"Le nombre d'angle n'est pas correct, il faut un total de {len(articulations)-1} angle(s).")
                return
        # Dans ces listes sont stockés les centres de masse
        Cx = [extr_socle["x"]]
        Cy = [extr_socle["y"]]
        # Pour toutes les articulations on calcule le centre de masse
        # (Pas 0 (socle) car déjà calculé dans les tableaux)
        for i in range(1,len(Ex)):
            temp_Cx = 0
            temp_Cy = 0
            # Boucle for pour faire la somme et calculer toutes les articulations
            for u in range(i+1):
                # Si on arrive à l'articulation où l'on cherche le centre de masse
                # on divise l'extrémité par 2 pour obtenir centre de l'articulation
                if i == u:
                    temp_Cx += Ex[u]/2
                    temp_Cy += Ey[u]/2
                    Cx.append(temp_Cx)
                    Cy.append(temp_Cy)
                    break
                # Sinon on additionne simplement les extrémités
                temp_Cx += Ex[u]
                temp_Cy += Ey[u]
        masse_totale = 0
        # Masse * postion pour tous les centres de masse
        pos_massX = 0
        pos_massY = 0
        for i in range(len(Cx)):
            masse_totale += articulations[str(i)]["masse"]
            pos_massX += articulations[str(i)]["masse"] * Cx[i]
            pos_massY += articulations[str(i)]["masse"] * Cy[i]
        # Applique la formule du centre de masse et arrondis à 4 chiffres après la virgule
        CM = (round(pos_massX/masse_totale, 4), round(pos_massY/masse_totale, 4))
        logger.success(f"Le centre de masse de la grue trouvé est {CM}[m]")
        return CM
    except KeyError as e:
        logger.error(f"Vous avez fait une erreur dans le fichier toml, Python le renseigne l'index {e}")
        return



def hauteur_enfoncement(masses, longueur_barge):
    """ 
    Calcule la hauteur de l'enfoncement de la plateforme
    Pré:    masses -> list
            longueur_barge -> int
    Post: retourne un int (hauteur de l'enfoncement)
    """
    somme_masses = sum(masses)
    return somme_masses / (DENSITE_EAU * longueur_barge**2)

def angles_immersion_soulevement(hauteur_enfoncement, hauteur_barge, longueur_barge):
    """
    Calcule les angles d'immersion et de soulèvement
    Pré:    longueur_barge -> int
            hauteur_enfoncement -> int
            hauteur_barge -> int
    Post: retourne une liste avec les deux angles
    """
    angle_immersion = math.atan((2*(hauteur_barge-hauteur_enfoncement))/longueur_barge)
    angle_soulevement = math.atan((2*hauteur_enfoncement)/longueur_barge)
    return [angle_immersion, angle_soulevement]

def aire_immergee(longueur_barge, hauteur_enfoncement):
    """
    Calcule simplement l'aire immergée sous l'eau
    Pré:    longueur_barge -> int
            hauteur_enfoncement -> int
    Post: retourne un int (aire)
    """
    return longueur_barge * hauteur_enfoncement

def longueurs_immergees(aire_immergee, longueur_barge, angles):
    """
    Calcule la longueur du côté gauche et droit qui est immergée
    Pré:    aire_immergee -> int
            longueur_barge -> int
            angles -> list
    Post: retourne une liste: [longueur_gauche, longueur_droite]
    """
    min_angle = min(angles)
    hauteur_gauche = aire_immergee/longueur_barge + (longueur_barge*math.tan(min_angle))/2
    hauteur_droite = aire_immergee/longueur_barge - (longueur_barge*math.tan(min_angle))/2
    return [hauteur_gauche, hauteur_droite]




centre_masse_grue(data["Grue"], [math.pi/6,0, -math.pi/6,0])

xpoints = np.array([0, 7])
ypoints = np.array([3, 3])
plt.plot(xpoints, ypoints, '--')
#plt.show()