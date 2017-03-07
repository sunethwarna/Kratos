// Author: Salva Latorre, latorre@cimne.upc.edu

#ifndef PERFORMANCE_UTILITIES_H
#define PERFORMANCE_UTILITIES_H

namespace Kratos {
    
    class PerformanceUtilities {

    public:

        typedef ModelPart::ElementsContainerType ElementsArrayType;
        typedef ModelPart::NodesContainerType NodesContainerType;

        PerformanceUtilities() {};

        virtual ~PerformanceUtilities() {};

        void PerformanceAccessingGetYoung(ModelPart& r_model_part) {
            double a = 0.0;
            for (ModelPart::ElementIterator i_elem = r_model_part.ElementsBegin(); i_elem != r_model_part.ElementsEnd(); ++i_elem) {
                Kratos::SphericParticle& particle = static_cast<Kratos::SphericParticle&>(*i_elem);
                const double young_modulus = particle.GetYoung();
                const double poisson = particle.GetPoisson();
                a += young_modulus + poisson;
            }
            if (a > 1.0) return;
        } //PerformanceAccessingToGetYoung
        
        void PerformanceAccessingGetPropertiesYoungModulus(ModelPart& r_model_part) {
            double a = 0.0;
            for (ModelPart::ElementIterator i_elem = r_model_part.ElementsBegin(); i_elem != r_model_part.ElementsEnd(); ++i_elem) {
                Kratos::SphericParticle& particle = static_cast<Kratos::SphericParticle&>(*i_elem);
                const double young_modulus = particle.GetProperties()[YOUNG_MODULUS];
                const double poisson = particle.GetProperties()[POISSON_RATIO];
                a += young_modulus + poisson;
            }
            if (a > 1.0) return;
        } //PerformanceAccessingToGetPropertiesYoungModulus
        
        void PerformanceAccessingProcessInfo(ModelPart& r_model_part) {

            for (ModelPart::ElementIterator i_elem = r_model_part.ElementsBegin(); i_elem != r_model_part.ElementsEnd(); ++i_elem) {
                Kratos::SphericParticle& particle = static_cast<Kratos::SphericParticle&>(*i_elem);
                const double time = r_model_part.GetProcessInfo()[TIME];
            }
        } //PerformanceAccessingToProcessInfo

    }; // Class PerformanceUtilities

} // namespace Kratos

#endif // PERFORMANCE_UTILITIES_H
