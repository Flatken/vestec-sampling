# Checkliste für die Anwendungsklasse 1 (Markdown-Template)
> Diese Checkliste liegt in der Version 1.0 vor. Sie basiert auf dem Dokument [QMH-DLR-04-V03-Anhang in Version 1.0](http://portal.dlr.de/Qualitaetsmanagement/QMH/QMH/Teil2/QMH-DLR-04-V03-Anhang.pdf).

## Nutzungshinweise
Diese Checkliste liefert Empfehlungen zur Software-Entwicklung. Primär richten sich diese an Software-Entwickler zur Selbsteinschätzung entwickelter Software und als Ideengeber für die weitere Entwicklung. Die Checkliste liefert keine neuen, revolutionären Ansätze zur Software-Entwicklung. Sie hilft aber notwendige, wesentliche Schritte der Software-Entwicklung nicht zu vergessen. Zudem können die Empfehlungen als Argumentationshilfe dienen.

Die Empfehlungen sind mit Fokus auf Wissenserhalt und gute Software-Engineering-Praxis erstellt. Sie unterstützen dabei, die Nachvollziehbarkeit und Nachhaltigkeit entwickelter Software zu erhalten. Die Empfehlungen motivieren den Einsatz von Tools, die Erstellung von Dokumentation, die Etablierung von Prozessen oder die Einhaltung von Prinzipien. Bei der Bewertung einer Empfehlung empfiehlt es sich daher zu überlegen, inwieweit der genannte Aspekt umgesetzt ist und ob Verbesserungsbedarf besteht. Dies kann man beispielsweise wie folgt umsetzen:

* Gibt es derzeit keinen Verbesserungsbedarf und die Empfehlung ist prinzipiell passend adressiert? Status: **ok**
* Gibt es Verbesserungspotential, welches in nächster Zeit umgesetzt werden kann bzw. sollte? Status: **todo**, Verbesserungsbedarf unter Bemerkungen festhalten
* Ist die Empfehlung derzeit noch nicht relevant, könnte aber in einer späteren Entwicklungsphase hilfreich sein? Status: **future**
* Ist die Empfehlung im Entwicklungskontext nicht sinnvoll umsetzbar? Status: **n.a.** (not applicable, nicht zutreffend), Grund unter Bemerkungen erläutern

> Den Status zwischen `**[]**` vermerken. Die Bemerkungen unterhalb der Empfehlung als Liste (z.B. `* Repository einrichten`) erfassen.

Zusätzliche Erläuterungen zu den Empfehlungen finden Sie im Dokument QMH-DLR-04-V03-Anhang. Zudem können Sie sich bei [Software-Engineering-Ansprechpartner ihres Instituts bzw. ihrer Einrichtung](http://portal.dlr.de/Organisation/Organisationseinheiten/ProgrDirektionQuerschnittsaufg/IKM/SWEngineering/Lists/Ansprechpersonen/AllItems.aspx) wenden.

## Inhaltsverzeichnis
[[Qualifizierung](#qualifizierung)] [[Anforderungsmanagement](#anforderungsmanagement)] [[Software-Architektur](#software-architektur)] [[Änderungsmanagement](#aenderungsmanagement)] [[Design und Implementierung](#design-implementierung)] [[Software-Test](#software-test)] [[Release-Management](#release-management)] [[Automatisierung und Abhängigkeitsmanagement](#automatisierung-abhaengigkeiten)] 

## Qualifizierung <a name="qualifizierung"></a>
**[ok]** Der Software-Verantwortliche kennt die verschiedenen Anwendungsklassen und weiß, welche für seine Software anzustreben ist. *(EQA.1, ab Anwendungsklasse 1)*

**[ok]** Der Software-Verantwortliche weiß, wie er gezielt Unterstützung zu Beginn und im Verlauf der Entwicklung anfordern und sich mit anderen Kollegen zum Thema Software-Entwicklung austauschen kann. *(EQA.2, ab Anwendungsklasse 1)*

**[ok]** Die an der Entwicklung Beteiligten ermitteln den Qualifikationsbedarf in Bezug auf ihre Rolle  und die angestrebte Anwendungsklasse. Sie kommunizieren diesen Bedarf an den Vorgesetzten. *(EQA.3, ab Anwendungsklasse 1)*

**[ok]** Den an der Entwicklung Beteiligten stehen die für ihre Aufgaben benötigten Werkzeuge zur Verfügung und sie sind geschult in deren Benutzung. *(EQA.4, ab Anwendungsklasse 1)*

## Anforderungsmanagement <a name="anforderungsmanagement"></a>
**[ok]** Die Aufgabenstellung ist mit allen Beteiligten abgestimmt und dokumentiert. Sie beschreibt in knapper, verständlicher Form die Ziele, den Zweck der Software, die wesentlichen Anforderungen und die angestrebte Anwendungsklasse. *(EAM.1, ab Anwendungsklasse 1)*

**[todo]** Die Randbedingungen sind erfasst. *(EAM.3, ab Anwendungsklasse 1)*

## Software-Architektur <a name="software-architektur"></a>
**[todo]** Wesentliche Architekturkonzepte und damit zusammenhängende Entscheidungen sind zumindest in knapper Form dokumentiert. *(ESA.2, ab Anwendungsklasse 1)*

## Änderungsmanagement <a name="aenderungsmanagement"></a>
**[todo]** Die wichtigsten Informationen, um zur Entwicklung beitragen zu können, sind an einer zentralen Stelle abgelegt. *(EÄM.2, ab Anwendungsklasse 1)*

**[ok]** Bekannte Fehler, wichtige ausstehende Aufgaben und Ideen sind zumindest stichpunktartig in einer Liste festgehalten und zentral abgelegt. *(EÄM.5, ab Anwendungsklasse 1)*

**[ok]** Ein Repository ist in einem Versionskontrollsystem eingerichtet. Das Repository ist angemessen strukturiert und enthält möglichst alle Artefakte, die zum Erstellen einer nutzbaren Version der Software und deren Test erforderlich sind. *(EÄM.7, ab Anwendungsklasse 1)*

**[ok]** Jede Änderung des Repository dient möglichst einem spezifischen Zweck, enthält eine verständliche Beschreibung und hinterlässt die Software möglichst in einem konsistenten, funktionierenden Zustand. *(EÄM.8, ab Anwendungsklasse 1)*

## Design und Implementierung <a name="design-implementierung"></a>
**[ok]** Es werden die üblichen Konstrukte und Lösungsansätze der gewählten Programmiersprache eingesetzt sowie ein Regelsatz hinsichtlich des Programmierstils konsequent angewendet. Der Regelsatz bezieht sich zumindest auf die Formatierung und Kommentierung. *(EDI.1, ab Anwendungsklasse 1)*

**[ok]** Die Software ist möglichst modular strukturiert. Die Module sind lose gekoppelt, d.h., ein einzelnes Modul hängt möglichst gering von anderen Modulen ab. *(EDI.2, ab Anwendungsklasse 1)*

**[ok]** Im Quelltext und in den Kommentaren sind möglichst wenig duplizierte Informationen enthalten. („Don`t repeat yourself.“) *(EDI.9, ab Anwendungsklasse 1)*

**[ok]** Es werden einfache, verständliche Lösungen bevorzugt eingesetzt.  („Keep it simple and stupid.“). *(EDI.10, ab Anwendungsklasse 1)*

## Software-Test <a name="software-test"></a>
**[todo]** Die grundlegenden Funktionen und Eigenschaften der Software werden in einer möglichst betriebsnahen Umgebung getestet. *(EST.4, ab Anwendungsklasse 1)*

**[ok]** Das Repository enthält möglichst alle für den Test der Software erforderlichen Artefakte. *(EST.10, ab Anwendungsklasse 1)*

## Release-Management <a name="release-management"></a>
**[ok]** Jedes Release besitzt eine eindeutige Release-Nummer. Anhand der Release-Nummer lässt sich der zugrunde liegende Softwarestand im Repository ermitteln. *(ERM.1, ab Anwendungsklasse 1)*

**[todo]** Das Release-Paket enthält oder verweist auf die Nutzer-Dokumentation. Sie besteht zumindest aus Installations-, Nutzungs- und Kontaktinformationen sowie den Release Notes. Im Fall der Weitergabe des Release-Pakets an Dritte außerhalb des DLR, sind die Lizenzbedingungen unbedingt beizulegen. *(ERM.2, ab Anwendungsklasse 1)*

**[ok]** Während der Release-Durchführung werden alle vorgesehenen Testaktivitäten ausgeführt. *(ERM.6, ab Anwendungsklasse 1)*

**[future]** Vor der Weitergabe des Release-Pakets an Dritte außerhalb des DLR ist sicherzustellen, dass eine Lizenz festgelegt ist, die Lizenzbestimmungen verwendeter Fremdsoftware eingehalten werden und alle erforderlichen Lizenzinformationen dem Release-Paket beigelegt sind. *(ERM.9, ab Anwendungsklasse 1)*

**[future]** Vor der Weitergabe des Release-Pakets an Dritte außerhalb des DLR ist sicherzustellen, dass die Regelungen zur Exportkontrolle eingehalten werden. *(ERM.10, ab Anwendungsklasse 1)*

## Automatisierung und Abhängigkeitsmanagement <a name="automatisierung-abhaengigkeiten"></a>
**[ok]** Der einfache Build-Prozess läuft grundlegend automatisiert ab und notwendige manuelle Schritte sind beschrieben. Zudem sind ausreichend Informationen zur Betriebs- und Entwicklungsumgebung vorhanden. *(EAA.1, ab Anwendungsklasse 1)*

**[ok]** Die Abhängigkeiten zum Erstellen der Software sind zumindest mit dem Namen, der Versionsnummer, dem Zweck, den Lizenzbestimmungen und der Bezugsquelle beschrieben. *(EAA.2, ab Anwendungsklasse 1)*

**[ok]** Das Repository enthält möglichst alle Bestandteile, um den Build-Prozess durchführen zu können. *(EAA.10, ab Anwendungsklasse 1)*


> Diese Checkliste liegt in der Version 1.0 vor. Sie basiert auf dem Dokument [QMH-DLR-04-V03-Anhang in Version 1.0](http://portal.dlr.de/Qualitaetsmanagement/QMH/QMH/Teil2/QMH-DLR-04-V03-Anhang.pdf).
