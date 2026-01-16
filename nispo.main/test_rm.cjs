const fs = require('fs');
const path = require('path');
const dist = path.resolve('./dist');
console.log(`Trying to remove: ${dist}`);
try {
    if (fs.existsSync(dist)) {
        fs.rmSync(dist, { recursive: true, force: true });
        console.log('Successfully removed dist');
    } else {
        console.log('dist does not exist');
    }
} catch (e) {
    console.error('Failed to remove dist:', e);
}
