#ifndef ____useGoodiNames_h
#define ____useGoodiNames_h
namespace Goodi {

// A do nothing routine that "uses" the Goodi *ExtensionName C string
// arrays so that compiler warnings under gcc go away.

   inline void useGooidNames() {
      (void)(EventExtensionName);
      (void)(GtiExtensionName);
      (void)(SpectrumExtensionName);
      (void)(EboundsExtensionName);
      (void)(SCExtensionName);
      (void)(EAExtensionName);
      (void)(PSFExtensionName);
      (void)(REDISExtensionName);
      (void)(LCExtensionName);
      (void)(MatrixExtensionName);
   }
} // namespace Goodi

#endif // ____useGoodiNames_h
