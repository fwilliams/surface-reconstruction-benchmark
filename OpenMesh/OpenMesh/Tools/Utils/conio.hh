#ifndef OPENMESH_UTILS_CONIO_HH
#define OPENMESH_UTILS_CONIO_HH
// ----------------------------------------------------------------------------
namespace OpenMesh {
namespace Utils {
// ----------------------------------------------------------------------------

/** Check if characters a pending in stdin.
 *
 *  \return Number of characters available to read.
 *
 *  \see getch(), getche()
 */
int kbhit(void);


/** A blocking single character input from stdin
 *
 *  \return Character, or -1 if an input error occurs.
 *
 *  \see getche(), kbhit()
 */
int getch(void);

/** A blocking single character input from stdin with echo.
 *
 *  \return Character, or -1 if an input error occurs.
 *  \see getch(), kbhit()
 */
int getche(void);

// ----------------------------------------------------------------------------
} // namespace Utils
} // namespace OpenMesh
// ----------------------------------------------------------------------------
#endif // OPENMESH_UTILS_CONIO_HH
// ============================================================================
