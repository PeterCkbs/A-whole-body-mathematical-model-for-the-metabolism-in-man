# athre:
# - Figure 6 seems to not be eps/pdf
    #Thank you, now as pdf.
# - In figure 9 (a) the title says flux (which is a rate, like in fig. 9 (b)), but the unit on the y-axis is [g]. In the text you write cumulative, so I suppose it is a sum of fluxes (or integral if it is continuous) which means [g] is correct, but make sure that the title, the y-axis and figure text match. Typically, we don't use titles on the plots and just describe them in the caption. 
# - Figure 10 and 11, maybe also 8, seem to have small text. As a rule of thump the size of the text on the plots should be the same as the size of the normal text. 
# - Remember to treat equations as if they are a normal sentence where e.g. you use commas when you list multiple equations or a dot when the equation is the end of a sentence.

# tobk:
# General comments
# ----------------
# - make sure that there are no errors when compiling the project. Unfortunately, Overleaf produces a pdf even though there are errors. But there might still be errors in the pdf which can be hard to spot.
# - we also try to keep the Latex code as simple as possible to avoid issues in the publishing process. For instance, Table II looks very nice, but it might be an advantage to 1) create an .eps figure of it and 2) insert that .eps into the paper instead of having the code directly there.
# - we recommend students to structure their Latex code in multiple files and folders (as you have done with the diagrams) to avoid long .tex files. For instance, each section might be a .tex file on its own, but long sections could also be split into multiple .tex files.
# - Some passages contain a fair amount of pathos, e.g., "... the models beauty shines, ..." More neutral or plain language is preferred for scientific papers.
# - John probably already told you to write active sentences, e.g., "We model something ..." instead of "Something is modeled ..." Occasionally, passive sentences are alright for variation or in specific situations, e.g., "Accurate modeling of the entire human body is considered a difficult task".
#
# Rough layout of the introduction
# --------------------------
# First part:
# briefly answer what and why.
#
# Second part :
# Elaborate on central concepts necessary to understand the contribution of this work.
#
# Third part:
# Describe the context of this work, i.e., describe what others have done.
#
# Fourth part:
# Describe what we do in this contribution.
#
# Fifth part:
# Outline of the paper.
#
# First part: In a single paragraph, you need to convince the reader that this paper is worth reading using high-level arguments. Why is the metabolism important - what is the problem? How can mathematical models help address this problem. What's special about whole-body multi-scale models? (You are already in the right direction when you say that a lot of reactions are happening in the human body).
#
# Second part: Essentially, this part should explain the reader enough for them to understand the third and fourth part, and you can also elaborate on why the subject is important.
#
# Third part: In this part, you cite other peoples work. It's important to be selective, i.e., choose papers that are relevant for the present work. Otherwise, it may confuse the reader. Towards the end of this part, you can summarize what _hasn't_ been done before such that the reader is ready for the description of your contribution in the fourth part.
#
# Fourth part: This part should start with a concise description of this works contribution (essentially, what is new) and provide a relevant amount of details such that the reader understands what the overall work is.
#
# Fifth part: This part is short and precise. It simply summarizes the contents of the remaining sections of the paper. It often starts with a variation of the sentence "The remaining part of this paper is structured as follows."


%%%%%%%%%%First part
An extreme number of reactions constantly happen in the body, to accommodate the 

Modeling all reactions mathematically is a challenging task 


%%%%%%%%%% Third part
Complex mathematical models exist in today's literature \cite{dash_li_kim_saidel_cabrera_2008,panunzi_pompa_borri_piemonte_gaetano_2020, sorensen_1978, yasemi_jolicoeur_2021} to describe the metabolic reactions in man, although a simple and intuitive framework, in which these advanced mathematical equations can easily be incorporated is not readily available. This type of modelling is specifically useful in PK-PD drug development, such as enzyme inhibitors/activators or hormonal drugs. To model these reactions in a physiologically relevant setting, requires us to look at the system as a whole-body model. A key benefit to looking at the system of metabolic reactions as a whole-body model, is that it allows us to simulate important metabolic concentrations throughout different compartments (organs) under various conditions. The system being affected by various metabolic concentrations of the different organs, such that homeostasis is attempted for all organs, and not just a single compartment. \\



%%% udkast
An extreme number of reactions happens in the body, where the formation and breakdown of new metabolites constantly happen to accommodate the different stages that the human body can be in. This metabolism is diverse, and not two people have the same pharmacokinetical parameters, which poses a challenge in PKPD drug development. However, this issue can be met using complex mathematical equations and parameters that are fitted to a general population. Doing so, can provide an expected concentration of a single metabolite in a single organ. It is necessary to realize that not two organs are the same, and such we need to be organ specific in our reactions. To model these reactions in a physiologically relevant setting, requires us to look at the system as a whole-body model. A key benefit to looking at the system of metabolic reactions as a whole-body model, is that it allows us to simulate important metabolic concentrations throughout different compartments (organs) under various conditions. Complex mathematical models exist in today's literature \cite{dash_li_kim_saidel_cabrera_2008,panunzi_pompa_borri_piemonte_gaetano_2020, sorensen_1978, yasemi_jolicoeur_2021} that describes the metabolic reactions in man, although a simple and intuitive framework, in which these advanced mathematical equations can easily be incorporated is not readily available. 