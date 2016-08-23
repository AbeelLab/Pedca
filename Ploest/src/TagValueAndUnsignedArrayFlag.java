
public class TagValueAndUnsignedArrayFlag {
    public final Object value;
    public final boolean isUnsignedArray;

    public TagValueAndUnsignedArrayFlag(Object value, boolean unsignedArray) {
        this.value = value;
        isUnsignedArray = unsignedArray;
    }

    TagValueAndUnsignedArrayFlag(Object value) {
        this.value = value;
        this.isUnsignedArray = false;
    }
}
