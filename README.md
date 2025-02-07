# Arsenic-Rational-Sampling-Design-in-Iowa
Research Project;  Part of the Iowa Grants-to-Counties (GTC) program.

Project Status: Completed;

Introduction to the paper "Risk Based Arsenic Rational Sampling Design for Public and Environmental Health Management", [link of the paper](https://arxiv.org/abs/2102.11118)

which has been accepted by *Chemometrics and Intelligent Laboratory Systems*

Auther: Lihao Yin(the first auther), Huiyan Sang, Douglas J. Schnoebelenc, Brian Wels, Don Simmons,Alyssa Mattson, Michael Schueller, Michael Pentelladan, Susie Y. Dai;

key words: spatially clustered function model; Graphic fused lasso; ADMM algorithm;

## Abstract
Groundwater contaminated with arsenic has been recognized as a global threat, which negatively impacts human health. Populations that rely on private wells for their drinking water are vulnerable to the potential arsenic-related health risks such as cancer and birth defects. Arsenic exposure through drinking water is among one of the primary arsenic exposure routes that can be effectively managed by active testing and water treatment. From the public and environmental health management perspective, it is critical to allocate the limited resources to establish an effective arsenic sampling and testing plan for health risk mitigation. We present a spatially adaptive sampling design approach, based on a spatially varying estimation of the underlying contamination distribution. 
 The method is different from traditional sampling design methods that often reply on a spatially constant or smoothly varying contamination distribution. 
 In contrast, we propose a statistical regularization method to automatically detect spatial clusters of the underlying contamination risk from the currently available private well arsenic testing data in the USA, Iowa. This allows us to develop a sampling design method that is adaptive to the changes in the contamination risk across the identified regions. 
 We provide the spatially adaptive sample size calculation and sampling locations determination for different acceptance precision and confidence levels for each cluster, to effectively mitigate the arsenic risk from the resource management perspectives. The model presents a framework that can be widely used for other environmental contaminant monitoring and sampling for public and environmental health. 
 
 ## Brief Introduction to the Project
 ### Data Description
 The raw data amount to 14,570 previously collected observations of Arsenic tests in total, and the Figure below shows the spatial distribution of the observations. The red points indicate the locations of the wells which are not contaiminated by Arsenic, and the blue points indicate the wells which are contaminated.
 ![image](Figures/overview1.png)
 
 The Figure below shows the locations of more than 400,000 wells in Iowa, from which we will sample to assess the water quality in the future.
 ![image](Figures/candidate.png)
 
 ### Project Aims
 To cluster Iowa state into several sub-regions and the risk of Arsenic contamination is the same within each subregion and different between any two subregions;
 
 To estimate the Arsenic risk in each sub-regions after clustering;
 
 To develope an economic sampling strategy to monitor the Arsenic contamination risk in Iowa; 
 
 ### Methodology
 Proposed a logistic regression model with spatially varying parameters, to model the spatially inhomogeneous Arsenic contamination risk in Iowa;
 
 Implemented the regularized tech with graphic fused Lasso to fit the proposed model;
 
 Proposed a thinning algorithm for efficiently sampling;
 
 ### Results
 In the Figure below, we divided Iowa into 3 clusters by the proposed methodologies. The Arsenic risks are 0.0287, 0.2088 and 0.3373 in cluster 1, 2, 3 respectively.
 ![image](Figures/ES6.png)
 
 The Figure below shows the sampling locations for water assessment in Iowa;
 ![image](Figures/NEW1.png)
 
 
