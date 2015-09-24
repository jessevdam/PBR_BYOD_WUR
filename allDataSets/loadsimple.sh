#stardog-admin --server snarl://ssb5:8080 db drop main -u admin -p ssbdata 
stardog-admin --server snarl://localhost:8080 db create -u admin -p ssbdata -i -n byodplant2 -o strict.parsing=false -- plantfixed.n3 snps.n3 accession_ontology.ttl accession.ttl byod_plant_trait_data.owl byod_plant_traits.owl plantontology.owl


