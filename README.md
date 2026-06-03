# AmbientTomo – Synthetic Aperture Sonar & Tomographie

Das primäre Ziel dieses Projekts ist die Entwicklung eines bildgebenden Synthetic Aperture Sonar (SAS) unter Einsatz moderner tomographischer Rekonstruktionstechniken. Aktuell konzentriert sich die Entwicklung auf die Implementierung und Validierung der **Vorwärtsmodellierung** (Forward Modeling).

## Was passiert hier?

* **Vorwärtsmodellierung (Forward Modeling):** Um ein Sonarsystem oder tomographische Algorithmen zu testen, muss die physikalische Realität präzise im Computer simuliert werden. Das System berechnet hierzu die Ausbreitung akustischer Wellen (vorerst) in einem zweidimensionalen Medium. Ein künstliches Signal (Puls) regt das Medium an, breitet sich aus, bricht sich an Materialgrenzen und erzeugt Reflexionen. Diese synthetischen Daten bilden die essentielle Grundlage für die spätere Rückrechnung (Inversion) in ein hochauflösendes Bild.
* **Absorbierende Randbedingungen:** Um unendliche Ausdehnungen (wie den offenen Ozean oder tiefe geologische Strukturen) auf einem begrenzten digitalen Gitternetz fehlerfrei darzustellen, fangen spezialisierte Ränder die Wellenfronten ab. Künstliche Echos an den Rändern des Simulationsbereichs werden geschluckt, sodass die Signale sauber aus dem Bild laufen und die Sensordaten nicht verfälschen. Die Simulation nutzt die Finite-Differenzen-Methode (FDTD) 
auf einem versetzten Gitter (Staggered Grid) und implementiert hochwirksame 
absorbierende Randbedingungen (Convolutional Perfectly Matched Layer - CPML).
* **Visualisierung & Dynamikregelung:** Da die akustischen Amplituden mit zunehmender Entfernung von der Signalquelle stark abnehmen, sorgt eine automatische Verstärkungsregelung dafür, dass sowohl der initiale hochenergetische Puls als auch die extrem schwachen, tiefen Reflexionen in der visuellen Ausgabe exakt und kontrastreich dargestellt werden.

## Integration mit Bazel
Das gesamte System ist vollständig deterministisch aufgebaut. Jede Modifikation an der physikalischen Vorwärtsmodellierung oder den Parametern der Wellenausbreitung führt automatisch zu einer exakt reproduzierbaren Neuberechnung des Simulationsergebnisses. Unveränderte Stände werden ohne erneuten Rechenaufwand direkt aus dem globalen Build-Cache geladen.